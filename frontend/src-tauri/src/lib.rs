// Tauri backend for Ligora
// Spawns the Python backend process and bridges newline-delimited JSON IPC
// (stdin/stdout) to the frontend via `send_command` / `backend-event`.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;
use tauri::{AppHandle, Emitter, Manager};

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .invoke_handler(tauri::generate_handler![
            send_command,
            get_status,
            copy_file
        ])
        .setup(|app| {
            let handle = app.handle().clone();
            match start_backend(handle) {
                Ok(state) => {
                    app.manage(state);
                }
                Err(e) => {
                    eprintln!("Failed to start Python backend: {e}");
                    app.manage(BackendState::failed(format!(
                        "Failed to start Python backend: {e}"
                    )));
                }
            }
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}

// ----------------------------------------------------------------------
// IPC data types
// ----------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CommandMessage {
    pub command_id: String,
    #[serde(rename = "type")]
    pub command_type: String,
    #[serde(default)]
    pub payload: serde_json::Value,
    #[serde(default)]
    pub session_id: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CommandResponse {
    pub command_id: String,
    pub success: bool,
    #[serde(default)]
    pub data: Option<serde_json::Value>,
    #[serde(default)]
    pub error: Option<String>,
}

// ----------------------------------------------------------------------
// Backend process state
// ----------------------------------------------------------------------

#[derive(Clone)]
struct BackendState {
    sender: Arc<Mutex<Option<std::sync::mpsc::Sender<CommandMessage>>>>,
    pending:
        Arc<Mutex<HashMap<String, std::sync::mpsc::Sender<CommandResponse>>>>,
    error_reason: Arc<Mutex<Option<String>>>,
}

impl BackendState {
    fn failed(reason: String) -> Self {
        Self {
            sender: Arc::new(Mutex::new(None)),
            pending: Arc::new(Mutex::new(HashMap::new())),
            error_reason: Arc::new(Mutex::new(Some(reason))),
        }
    }
}

fn find_backend_script() -> Option<std::path::PathBuf> {
    // Explicit override for packaged/snap deployments.
    if let Ok(dir) = std::env::var("LIGORA_BACKEND_DIR") {
        let p = std::path::PathBuf::from(dir).join("server.py");
        if p.exists() {
            return Some(p);
        }
    }
    // Dev fallbacks: run straight from the repository checkout.
    let candidates = [
        std::path::PathBuf::from("../../backend/ligora_backend"),
        std::path::PathBuf::from("../../../backend/ligora_backend"),
        std::path::PathBuf::from("/usr/lib/ligora/backend"),
        std::path::PathBuf::from("/opt/ligora/backend"),
    ];
    for dir in candidates {
        let p = dir.join("server.py");
        if p.exists() {
            return Some(dir);
        }
    }
    None
}

fn start_backend(app: AppHandle) -> Result<BackendState, String> {
    let module_dir =
        find_backend_script().ok_or_else(|| {
            "Python backend (ligora_backend/server.py) not found".to_string()
        })?;

    let python =
        std::env::var("LIGORA_PYTHON").unwrap_or_else(|_| "python3".into());
    let module_dir_str = module_dir.to_string_lossy().to_string();

    let mut child = Command::new(python)
        .arg("-u") // unbuffered stdout for line-based IPC
        .arg("-m")
        .arg("ligora_backend.server")
        .current_dir(&module_dir)
        .env("PYTHONPATH", &module_dir_str)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("failed to spawn python backend: {e}"))?;

    let stdin = child
        .stdin
        .take()
        .ok_or_else(|| "failed to capture backend stdin".to_string())?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "failed to capture backend stdout".to_string())?;

    let (cmd_tx, cmd_rx) = std::sync::mpsc::channel::<CommandMessage>();
    let pending: Arc<Mutex<HashMap<String, std::sync::mpsc::Sender<CommandResponse>>>> =
        Arc::new(Mutex::new(HashMap::new()));

    // Writer: forwards queued commands to the backend's stdin.
    thread::spawn(move || {
        let mut writer = stdin;
        for message in cmd_rx {
            match serde_json::to_string(&message) {
                Ok(json) => {
                    if writer
                        .write_all(json.as_bytes())
                        .and_then(|_| writer.write_all(b"\n"))
                        .and_then(|_| writer.flush())
                        .is_err()
                    {
                        break;
                    }
                }
                Err(_) => break,
            }
        }
    });

    // Reader: correlates responses by command_id; forwards events.
    {
        let pending = pending.clone();
        let app_for_events = app.clone();
        thread::spawn(move || {
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                let line = match line {
                    Ok(l) => l,
                    Err(_) => break,
                };
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                if let Ok(response) =
                    serde_json::from_str::<CommandResponse>(trimmed)
                {
                    if let Some(tx) =
                        pending.lock().unwrap().remove(&response.command_id)
                    {
                        let _ = tx.send(response);
                    }
                    continue;
                }
                if let Ok(value) = serde_json::from_str::<serde_json::Value>(trimmed)
                {
                    if value.get("event").is_some() {
                        let _ = app_for_events.emit("backend-event", value);
                    }
                }
            }
        });
    }

    // Reaper: reap the child when it exits (avoids zombie processes).
    thread::spawn(move || {
        let _ = child.wait();
    });

    Ok(BackendState {
        sender: Arc::new(Mutex::new(Some(cmd_tx))),
        pending,
        error_reason: Arc::new(Mutex::new(None)),
    })
}

// ----------------------------------------------------------------------
// Tauri commands
// ----------------------------------------------------------------------

#[tauri::command]
async fn send_command(
    app: AppHandle,
    command_type: String,
    payload: serde_json::Value,
    session_id: String,
) -> Result<CommandResponse, String> {
    let state = app.state::<BackendState>().inner().clone();

    let sender = state
        .sender
        .lock()
        .unwrap()
        .clone()
        .ok_or_else(|| {
            state
                .error_reason
                .lock()
                .unwrap()
                .clone()
                .unwrap_or_else(|| "Python backend is not running".into())
        })?;

    let message = CommandMessage {
        command_id: uuid::Uuid::new_v4().to_string(),
        command_type,
        payload,
        session_id,
    };

    // Register the responder BEFORE writing so the reader can find it.
    let (tx, rx) = std::sync::mpsc::channel::<CommandResponse>();
    state
        .pending
        .lock()
        .unwrap()
        .insert(message.command_id.clone(), tx);

    sender
        .send(message.clone())
        .map_err(|e| format!("failed to send command: {e}"))?;

    // Long timeout: docking jobs legitimately run for minutes.
    match rx.recv_timeout(Duration::from_secs(3600)) {
        Ok(response) => Ok(response),
        Err(_) => {
            state.pending.lock().unwrap().remove(&message.command_id);
            Err("Python backend did not respond in time".to_string())
        }
    }
}

#[tauri::command]
async fn get_status(app: AppHandle) -> Result<serde_json::Value, String> {
    let resp = send_command(
        app,
        "get_status".into(),
        serde_json::json!({}),
        "".into(),
    )
    .await?;
    Ok(resp.data.unwrap_or_else(|| serde_json::json!({})))
}

#[tauri::command]
async fn copy_file(src: String, dest: String) -> Result<bool, String> {
    let data = std::fs::read(&src).map_err(|e| format!("read failed: {e}"))?;
    if let Some(parent) = std::path::Path::new(&dest).parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    std::fs::write(&dest, data).map_err(|e| format!("write failed: {e}"))?;
    Ok(true)
}
