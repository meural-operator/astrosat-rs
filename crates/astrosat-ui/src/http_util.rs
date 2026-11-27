//! HTTP utilities, query/URL parsing, safe path resolution, and JSON response helpers.

use std::collections::HashMap;
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use serde::Serialize;
use tiny_http::{Header, Request, Response, StatusCode};

/// Maximum JSON request body accepted by POST endpoints (4 MiB).
pub const MAX_BODY_BYTES: usize = 4 * 1024 * 1024;

/// Structured HTTP error that converts cleanly into a JSON `{ "error", "status" }` response.
#[derive(Debug)]
pub struct ApiError {
    pub status: u16,
    pub message: String,
}

impl ApiError {
    pub fn bad_request(msg: impl Into<String>) -> Self {
        Self {
            status: 400,
            message: msg.into(),
        }
    }

    pub fn forbidden(msg: impl Into<String>) -> Self {
        Self {
            status: 403,
            message: msg.into(),
        }
    }

    pub fn not_found(msg: impl Into<String>) -> Self {
        Self {
            status: 404,
            message: msg.into(),
        }
    }

    pub fn internal(msg: impl Into<String>) -> Self {
        Self {
            status: 500,
            message: msg.into(),
        }
    }

    pub fn method_not_allowed(msg: impl Into<String>) -> Self {
        Self {
            status: 405,
            message: msg.into(),
        }
    }
}

impl std::fmt::Display for ApiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HTTP {}: {}", self.status, self.message)
    }
}

impl std::error::Error for ApiError {}

pub type ApiResult<T> = std::result::Result<T, ApiError>;

#[derive(Serialize)]
struct ErrorBody<'a> {
    error: &'a str,
    status: u16,
}

fn header(name: &str, value: &str) -> Header {
    Header::from_bytes(name.as_bytes(), value.as_bytes()).expect("valid static ASCII HTTP header")
}

/// Build an HTTP 200 JSON response.
pub fn json_ok<T: Serialize>(value: &T) -> Response<std::io::Cursor<Vec<u8>>> {
    json_status(200, value)
}

/// Build a JSON response with an arbitrary HTTP status code.
pub fn json_status<T: Serialize>(status: u16, value: &T) -> Response<std::io::Cursor<Vec<u8>>> {
    match serde_json::to_vec(value) {
        Ok(bytes) => Response::from_data(bytes)
            .with_status_code(StatusCode(status))
            .with_header(header("Content-Type", "application/json; charset=utf-8"))
            .with_header(header("Cache-Control", "no-store"))
            .with_header(header("X-Content-Type-Options", "nosniff")),
        Err(e) => {
            let fallback = serde_json::to_vec(&ErrorBody {
                error: &format!("JSON serialization failed: {e}"),
                status: 500,
            })
            .unwrap_or_else(|_| b"{\"error\":\"serialization error\",\"status\":500}".to_vec());
            Response::from_data(fallback)
                .with_status_code(StatusCode(500))
                .with_header(header("Content-Type", "application/json; charset=utf-8"))
                .with_header(header("Cache-Control", "no-store"))
        }
    }
}

/// Build an HTTP error JSON response from an [`ApiError`].
pub fn json_error(err: &ApiError) -> Response<std::io::Cursor<Vec<u8>>> {
    json_status(
        err.status,
        &ErrorBody {
            error: &err.message,
            status: err.status,
        },
    )
}

/// Build a static asset response (HTML, CSS, JS).
pub fn static_response(
    content_type: &str,
    body: &'static str,
) -> Response<std::io::Cursor<Vec<u8>>> {
    Response::from_data(body.as_bytes().to_vec())
        .with_status_code(StatusCode(200))
        .with_header(header("Content-Type", content_type))
        .with_header(header("Cache-Control", "no-cache"))
        .with_header(header("X-Content-Type-Options", "nosniff"))
}

/// Split a request URL into `(decoded_path, query_params)`.
pub fn parse_url(raw_url: &str) -> (String, HashMap<String, String>) {
    let (path_raw, query_raw) = match raw_url.split_once('?') {
        Some((p, q)) => (p, Some(q)),
        None => (raw_url, None),
    };
    let path = percent_decode(path_raw);
    let mut params = HashMap::new();
    if let Some(q) = query_raw {
        for pair in q.split('&') {
            if pair.is_empty() {
                continue;
            }
            if let Some((k, v)) = pair.split_once('=') {
                params.insert(percent_decode_query(k), percent_decode_query(v));
            } else {
                params.insert(percent_decode_query(pair), String::new());
            }
        }
    }
    (path, params)
}

/// Percent-decode a URL path segment (does NOT turn `+` into space).
pub fn percent_decode(input: &str) -> String {
    percent_decode_inner(input, false)
}

/// Percent-decode a query parameter key or value (`+` -> space).
pub fn percent_decode_query(input: &str) -> String {
    percent_decode_inner(input, true)
}

fn percent_decode_inner(input: &str, plus_as_space: bool) -> String {
    let bytes = input.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' if plus_as_space => {
                out.push(b' ');
                i += 1;
            }
            b'%' if i + 2 < bytes.len() => {
                if let (Some(hi), Some(lo)) = (hex_val(bytes[i + 1]), hex_val(bytes[i + 2])) {
                    out.push((hi << 4) | lo);
                    i += 3;
                } else {
                    out.push(b'%');
                    i += 1;
                }
            }
            b => {
                out.push(b);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn hex_val(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

/// Read a bounded JSON request body and deserialize it into `T`.
pub fn read_json_body<T: serde::de::DeserializeOwned>(req: &mut Request) -> ApiResult<T> {
    let mut buf = Vec::new();
    let reader = req.as_reader();
    let mut limited = reader.take((MAX_BODY_BYTES + 1) as u64);
    limited
        .read_to_end(&mut buf)
        .map_err(|e| ApiError::bad_request(format!("failed to read request body: {e}")))?;
    if buf.len() > MAX_BODY_BYTES {
        return Err(ApiError::bad_request(format!(
            "request body exceeds {} bytes limit",
            MAX_BODY_BYTES
        )));
    }
    if buf.is_empty() {
        return Err(ApiError::bad_request("empty JSON request body"));
    }
    serde_json::from_slice(&buf)
        .map_err(|e| ApiError::bad_request(format!("invalid JSON body: {e}")))
}

/// Parse an optional boolean query parameter (`true`/`1`/`yes` vs `false`/`0`/`no`).
pub fn query_bool(params: &HashMap<String, String>, key: &str, default: bool) -> ApiResult<bool> {
    match params.get(key).map(|s| s.trim()) {
        None | Some("") => Ok(default),
        Some(v) if v.eq_ignore_ascii_case("true") || v == "1" || v.eq_ignore_ascii_case("yes") => {
            Ok(true)
        }
        Some(v) if v.eq_ignore_ascii_case("false") || v == "0" || v.eq_ignore_ascii_case("no") => {
            Ok(false)
        }
        Some(other) => Err(ApiError::bad_request(format!(
            "invalid boolean parameter '{key}': '{other}'"
        ))),
    }
}

/// Parse an optional `usize` query parameter.
pub fn query_usize(
    params: &HashMap<String, String>,
    key: &str,
    default: usize,
) -> ApiResult<usize> {
    match params.get(key).map(|s| s.trim()) {
        None | Some("") => Ok(default),
        Some(v) => v.parse::<usize>().map_err(|_| {
            ApiError::bad_request(format!("invalid integer parameter '{key}': '{v}'"))
        }),
    }
}

/// Parse a required non-empty query parameter.
pub fn query_required<'a>(params: &'a HashMap<String, String>, key: &str) -> ApiResult<&'a str> {
    match params.get(key).map(|s| s.trim()) {
        Some(v) if !v.is_empty() => Ok(v),
        _ => Err(ApiError::bad_request(format!(
            "missing required query parameter '{key}'"
        ))),
    }
}

/// Check if a raw user-supplied path contains traversal (`..`), NUL bytes, or sensitive system directories.
fn validate_raw_path(raw: &str) -> ApiResult<PathBuf> {
    if raw.is_empty() {
        return Err(ApiError::bad_request("empty path"));
    }
    if raw.contains('\0') {
        return Err(ApiError::bad_request("path contains NUL byte"));
    }
    // Reject Windows-style traversal or Unix traversal explicitly before canonicalization.
    let normalized_slashes = raw.replace('\\', "/");
    for seg in normalized_slashes.split('/') {
        if seg == ".." {
            return Err(ApiError::bad_request(
                "path traversal ('..') is not allowed",
            ));
        }
    }
    let p = PathBuf::from(raw);
    for comp in p.components() {
        if matches!(comp, Component::ParentDir) {
            return Err(ApiError::bad_request(
                "path traversal ('..') is not allowed",
            ));
        }
    }
    Ok(p)
}

fn is_forbidden_system_path(p: &Path) -> bool {
    let s = p.to_string_lossy();
    let lower = s.to_ascii_lowercase();
    const FORBIDDEN_PREFIXES: &[&str] = &[
        "/etc",
        "/proc",
        "/sys",
        "/dev",
        "/boot",
        "/root",
        "/var/run",
        "/run/secrets",
    ];
    for pref in FORBIDDEN_PREFIXES {
        if lower == *pref || lower.starts_with(&format!("{pref}/")) {
            return true;
        }
    }
    lower.contains("/.ssh") || lower.contains("\\.ssh")
}

/// Resolve a user-supplied file path safely (relative to `data_dir` if relative).
/// Rejects `..` traversal, NUL bytes, sensitive system directories, and missing files.
pub fn resolve_existing_file(raw: &str, data_dir: &Path) -> ApiResult<PathBuf> {
    let p = validate_raw_path(raw)?;
    let candidate = if p.is_absolute() {
        p
    } else {
        data_dir.join(p)
    };
    if is_forbidden_system_path(&candidate) {
        return Err(ApiError::forbidden(format!(
            "access denied for system path: {}",
            candidate.display()
        )));
    }
    let canon = std::fs::canonicalize(&candidate).map_err(|e| match e.kind() {
        std::io::ErrorKind::NotFound => {
            ApiError::not_found(format!("file not found: {}", candidate.display()))
        }
        std::io::ErrorKind::PermissionDenied => {
            ApiError::forbidden(format!("permission denied: {}", candidate.display()))
        }
        _ => ApiError::bad_request(format!("cannot access {}: {e}", candidate.display())),
    })?;
    if is_forbidden_system_path(&canon) {
        return Err(ApiError::forbidden(format!(
            "access denied for system path: {}",
            canon.display()
        )));
    }
    if !canon.is_file() {
        return Err(ApiError::bad_request(format!(
            "path is not a regular file: {}",
            canon.display()
        )));
    }
    Ok(canon)
}

/// Resolve a user-supplied directory path safely (empty -> `data_dir`).
pub fn resolve_existing_dir(raw: Option<&str>, data_dir: &Path) -> ApiResult<PathBuf> {
    let candidate = match raw.map(str::trim) {
        None | Some("") => data_dir.to_path_buf(),
        Some(s) => {
            let p = validate_raw_path(s)?;
            if p.is_absolute() {
                p
            } else {
                data_dir.join(p)
            }
        }
    };
    if is_forbidden_system_path(&candidate) {
        return Err(ApiError::forbidden(format!(
            "access denied for system directory: {}",
            candidate.display()
        )));
    }
    let canon = std::fs::canonicalize(&candidate).map_err(|e| match e.kind() {
        std::io::ErrorKind::NotFound => {
            ApiError::not_found(format!("directory not found: {}", candidate.display()))
        }
        std::io::ErrorKind::PermissionDenied => {
            ApiError::forbidden(format!("permission denied: {}", candidate.display()))
        }
        _ => ApiError::bad_request(format!("cannot access {}: {e}", candidate.display())),
    })?;
    if is_forbidden_system_path(&canon) {
        return Err(ApiError::forbidden(format!(
            "access denied for system directory: {}",
            canon.display()
        )));
    }
    if !canon.is_dir() {
        return Err(ApiError::bad_request(format!(
            "path is not a directory: {}",
            canon.display()
        )));
    }
    Ok(canon)
}

/// Resolve an output file or directory path safely (does not need to exist yet).
pub fn resolve_output_path(raw: &str, data_dir: &Path) -> ApiResult<PathBuf> {
    let p = validate_raw_path(raw)?;
    let candidate = if p.is_absolute() {
        p
    } else {
        data_dir.join(p)
    };
    if is_forbidden_system_path(&candidate) {
        return Err(ApiError::forbidden(format!(
            "refusing to write to system path: {}",
            candidate.display()
        )));
    }
    let s = candidate.to_string_lossy().to_ascii_lowercase();
    for pref in ["/usr", "/bin", "/sbin", "/lib", "/lib64"] {
        if s == pref || s.starts_with(&format!("{pref}/")) {
            return Err(ApiError::forbidden(format!(
                "refusing to write to system path: {}",
                candidate.display()
            )));
        }
    }
    Ok(candidate)
}
