from __future__ import annotations

import uvicorn


if __name__ == "__main__":
    uvicorn.run("brain_app.app:app", host="127.0.0.1", port=8001, reload=False)
