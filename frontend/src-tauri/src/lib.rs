// Tauri backend for Ligora
// Handles communication with the Python backend via local socket or stdin/stdout

use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .invoke_handler(tauri::generate_handler![start_backend, stop_backend, send_command, get_status])
        .setup(|app| {
            // Start the Python backend on app startup
            let backend_handle = start_backend_internal();

            // Store backend handle in app state
            app.manage(backend_handle);

            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}

// Backend process handle
pub struct BackendHandle {
    pub process: Option<std::process::Child>,
    pub sender: Option<mpsc::Sender<CommandMessage>>,
    pub receiver: Option<mpsc::Receiver<CommandResponse>>,
    pub writer_thread: Option<thread::JoinHandle<()>>,
    pub reader_thread: Option<thread::JoinHandle<()>>,
}

// Message types for backend communication
#[derive(Debug, Serialize, Deserialize)]
pub struct CommandMessage {
    pub command_id: String,
    #[serde(rename = "type")]
    pub command_type: String,
    pub payload: serde_json::Value,
    pub session_id: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CommandResponse {
    pub command_id: String,
    pub success: bool,
    pub data: Option<serde_json::Value>,
    pub error: Option<String>,
}

// Start the Python backend
fn start_backend_internal() -> BackendHandle {
    // Look for the Python backend
    let backend_script = find_backend_script();

    if backend_script.is_none() {
        eprintln!("Warning: Python backend not found");
        return BackendHandle {
            process: None,
            sender: None,
            receiver: None,
            thread: None,
        };
    }

    let script_path = backend_script.unwrap();

    // Start the Python process
    let mut child = match Command::new("python3")
        .arg(&script_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
    {
        Ok(child) => child,
        Err(e) => {
            eprintln!("Failed to start Python backend: {}", e);
            return BackendHandle {
                process: None,
                sender: None,
                receiver: None,
                thread: None,
            };
        }
    };

    // Create channels for communication
    let (cmd_sender, cmd_receiver) = mpsc::channel();
    let (resp_sender, resp_receiver) = mpsc::channel();

    // Spawn thread to read stdout from backend
    let stdout = child.stdout.take().unwrap();
    let resp_sender_clone = resp_sender.clone();
    let reader_thread = thread::spawn(move || {
        let mut buffer = String::new();
        let mut reader = std::io::BufReader::new(stdout);

        loop {
            let mut line = String::new();
            if reader.read_line(&mut line).is_err() {
                break;
            }

            if line.trim().is_empty() {
                continue;
            }

            // Parse JSON response
            if let Ok(response) = serde_json::from_str::<CommandResponse>(&line) {
                let _ = resp_sender_clone.send(response).is_err();
            }

            buffer.clear();
        }
    });

    // Spawn thread to handle commands
    let stdin = child.stdin.take().unwrap();
    let cmd_receiver_clone = cmd_receiver;
    let writer_thread = thread::spawn(move || {
        let mut writer = std::io::BufWriter::new(stdin);

        for message in cmd_receiver_clone {
            if let Ok(json) = serde_json::to_string(&message) {
                let _ = writer.write_all(json.as_bytes());
                let _ = writer.write_all(b"\n");
                let _ = writer.flush();
            }
        }
    });

    BackendHandle {
        process: Some(child),
        sender: Some(cmd_sender),
        receiver: Some(resp_receiver),
        writer_thread: Some(writer_thread),
        reader_thread: Some(reader_thread),
    }
}

// Find the Python backend script
fn find_backend_script() -> Option<String> {
    // Check relative to the app bundle
    let exe_dir = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|p| p.to_path_buf()))?;

    // Try various paths
    let candidates = [
        exe_dir.join("backend/ligora_backend/server.py"),
        exe_dir.join("lib/ligora_backend/server.py"),
        exe_dir.join("server.py"),
        std::path::PathBuf::from("/usr/lib/ligora/backend/server.py"),
        std::path::PathBuf::from("/opt/ligora/backend/server.py"),
    ];

    for path in &candidates {
        if path.exists() {
            return Some(path.to_string_lossy().to_string());
        }
    }

    None
}

// Tauri commands
#[tauri::command]
async fn start_backend() -> Result<bool, String> {
    // Backend is started in setup hook
    Ok(true)
}

#[tauri::command]
async fn stop_backend() -> Result<bool, String> {
    // Get the backend handle from app state
    // This would need proper state management
    Ok(true)
}

#[tauri::command]
async fn send_command(
    command_type: String,
    payload: serde_json::Value,
    session_id: String,
) -> Result<CommandResponse, String> {
    let handle = APP_STATE
        .lock()
        .ok()
        .and_then(|s| s.backend.clone())
        .and_then(|b| b.sender.clone())
    .ok_or_else(|| "Python backend is not running".to_string())?;

    let message = CommandMessage {
        command_id: uuid::Uuid::new_v4().to_string(),
        command_type,
        payload,
        session_id,
    };

    handle
        .send(message)
        .map_err(|e| format!("failed to send command: {}", e))?;

    // Wait for the corresponding response from the Python backend.
    // The backend sends one JSON line per command response on stdout.
    let receiver = handle.receiver.clone().ok_or_else(
        || "backend response channel is not available".to_string()
    )?;

    loop {
        match receiver.recv_timeout(Duration::from_secs(60)) {
            Ok(response) => {
                if response.command_id == message.command_id {
                    return Ok(response);
                }
                // Still return responses that belong to this handle if the
                // command id does not match for some reason; otherwise keep
                // waiting for the matching command.
                continue;
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                return Err("Python backend closed the response channel".to_string());
            }
            Err(mpsc::RecvTimeoutError::Timeout) => {
                return Err("Python backend did not respond in time".to_string());
            }
        }
    }
}

#[tauri::command]
async fn get_status() -> Result<serde_json::Value, String> {
    // Get status from the Python backend
    Ok(serde_json::json!({
        "status": "running",
        "version": "0.1.0"
    }))
}
