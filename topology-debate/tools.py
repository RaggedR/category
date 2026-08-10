"""Tools available to debate agents: Python execution, file read/write."""

import subprocess
import sys
from pathlib import Path

# Directories agents can access
PUZZLES_DIR = Path(__file__).resolve().parent.parent.parent / "puzzles"
WORK_DIR = Path(__file__).resolve().parent

# Allowed read paths (agents can read puzzles code + their own workspace)
ALLOWED_READ_ROOTS = [PUZZLES_DIR, WORK_DIR]

# Allowed write paths (agents can only write to their workspace)
ALLOWED_WRITE_ROOTS = [WORK_DIR]


def run_python(code: str, timeout: int = 30) -> dict:
    """Execute Python code in a subprocess. Returns stdout and stderr."""
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(PUZZLES_DIR),  # run from puzzles dir so imports work
        )
        return {
            "stdout": result.stdout[:4000],
            "stderr": result.stderr[:2000],
            "returncode": result.returncode,
        }
    except subprocess.TimeoutExpired:
        return {
            "stdout": "",
            "stderr": f"Execution timed out after {timeout}s",
            "returncode": -1,
        }


def read_file(path: str) -> dict:
    """Read a file. Path must be under puzzles/ or topology-debate/."""
    p = Path(path).resolve()
    if not any(p.is_relative_to(root) for root in ALLOWED_READ_ROOTS):
        return {"error": f"Access denied: {path} is outside allowed directories"}
    if not p.exists():
        return {"error": f"File not found: {path}"}
    try:
        content = p.read_text()
        if len(content) > 8000:
            content = content[:8000] + f"\n... [truncated, {len(content)} chars total]"
        return {"content": content}
    except Exception as e:
        return {"error": str(e)}


def write_file(path: str, content: str) -> dict:
    """Write a file. Path must be under topology-debate/."""
    p = Path(path).resolve()
    if not any(p.is_relative_to(root) for root in ALLOWED_WRITE_ROOTS):
        return {"error": f"Access denied: can only write to {WORK_DIR}"}
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
        return {"success": True, "path": str(p)}
    except Exception as e:
        return {"error": str(e)}


# Tool definitions for the Anthropic SDK
TOOL_DEFINITIONS = [
    {
        "name": "run_python",
        "description": (
            "Execute Python code. The working directory is ~/git/puzzles/ so you "
            "can import lr_coefficients, solve_puzzle, transfer_operators directly. "
            "Use this to test conjectures, enumerate puzzles, compute examples, "
            "and verify mathematical claims. Print results to stdout."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "code": {
                    "type": "string",
                    "description": "Python code to execute",
                },
            },
            "required": ["code"],
        },
    },
    {
        "name": "read_file",
        "description": (
            "Read a file from ~/git/puzzles/ or ~/git/category/topology-debate/. "
            "Use this to inspect source code, read existing results, or check "
            "puzzle data files."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Absolute path or path relative to ~/git/puzzles/",
                },
            },
            "required": ["path"],
        },
    },
    {
        "name": "write_file",
        "description": (
            "Write a file to ~/git/category/topology-debate/. Use this to save "
            "results, write analysis scripts, or record findings."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "Absolute path (must be under topology-debate/)",
                },
                "content": {
                    "type": "string",
                    "description": "File content to write",
                },
            },
            "required": ["path", "content"],
        },
    },
]


def execute_tool(name: str, input: dict) -> str:
    """Dispatch a tool call and return the result as a string."""
    import json

    if name == "run_python":
        result = run_python(input["code"])
    elif name == "read_file":
        result = read_file(input["path"])
    elif name == "write_file":
        result = write_file(input["path"], input["content"])
    else:
        result = {"error": f"Unknown tool: {name}"}

    return json.dumps(result, indent=2)
