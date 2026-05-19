from __future__ import annotations

from pathlib import Path
import sys


def main() -> None:
    try:
        from streamlit.web import cli as streamlit_cli
    except ImportError as exc:
        raise SystemExit(
            "Streamlit is not installed. Install the app extra with "
            '`python -m pip install -e ".[app]"`. For local ND2 files, use '
            '`python -m pip install -e ".[app-nd2]"`.'
        ) from exc

    app_path = Path(__file__).with_name("app.py")
    sys.argv = ["streamlit", "run", str(app_path), *sys.argv[1:]]
    raise SystemExit(streamlit_cli.main())


if __name__ == "__main__":
    main()
