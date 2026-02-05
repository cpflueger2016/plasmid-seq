#!/usr/bin/env python3
"""
Validate a plasmid sample report JSON against a JSON schema.

Usage:
  scripts/validate_sample_report.py \
    -i references/sample_report.example.json \
    -s references/sample_report.schema.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate sample report JSON.")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to sample_report.json to validate.",
    )
    parser.add_argument(
        "-s",
        "--schema",
        default="references/sample_report.schema.json",
        help="Path to JSON schema (default: references/sample_report.schema.json).",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_with_jsonschema(instance: dict, schema: dict) -> tuple[bool, str]:
    try:
        import jsonschema  # type: ignore
    except Exception:
        return validate_fallback(instance, schema)

    try:
        jsonschema.validate(instance=instance, schema=schema)
    except jsonschema.exceptions.ValidationError as exc:
        path = ".".join(str(p) for p in exc.absolute_path) or "<root>"
        return False, f"{path}: {exc.message}"
    except Exception as exc:
        return False, str(exc)
    return True, "valid"


def _matches_type(value: object, type_name: str) -> bool:
    if type_name == "object":
        return isinstance(value, dict)
    if type_name == "array":
        return isinstance(value, list)
    if type_name == "string":
        return isinstance(value, str)
    if type_name == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if type_name == "number":
        return (isinstance(value, int) or isinstance(value, float)) and not isinstance(value, bool)
    if type_name == "boolean":
        return isinstance(value, bool)
    return True


def _validate_node(value: object, schema: dict, path: str, errors: list[str]) -> None:
    type_name = schema.get("type")
    if isinstance(type_name, str) and not _matches_type(value, type_name):
        errors.append(f"{path}: expected type '{type_name}'")
        return

    if "enum" in schema and value not in schema["enum"]:
        errors.append(f"{path}: value not in enum {schema['enum']}")

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            errors.append(f"{path}: value {value} < minimum {schema['minimum']}")
        if "maximum" in schema and value > schema["maximum"]:
            errors.append(f"{path}: value {value} > maximum {schema['maximum']}")

    if isinstance(value, dict):
        required = schema.get("required", [])
        for key in required:
            if key not in value:
                errors.append(f"{path}.{key}: missing required key")
        props = schema.get("properties", {})
        for key, child in props.items():
            if key in value and isinstance(child, dict):
                _validate_node(value[key], child, f"{path}.{key}", errors)

    if isinstance(value, list):
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for idx, item in enumerate(value):
                _validate_node(item, item_schema, f"{path}[{idx}]", errors)


def validate_fallback(instance: dict, schema: dict) -> tuple[bool, str]:
    errors: list[str] = []
    _validate_node(instance, schema, "<root>", errors)
    if errors:
        return False, errors[0]
    return True, "valid (fallback validator)"


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    schema_path = Path(args.schema).expanduser()

    if not input_path.exists():
        print(f"[validate][ERROR] input not found: {input_path}", file=sys.stderr)
        return 1
    if not schema_path.exists():
        print(f"[validate][ERROR] schema not found: {schema_path}", file=sys.stderr)
        return 1

    try:
        instance = read_json(input_path)
    except Exception as exc:
        print(f"[validate][ERROR] failed to parse input JSON: {exc}", file=sys.stderr)
        return 1

    try:
        schema = read_json(schema_path)
    except Exception as exc:
        print(f"[validate][ERROR] failed to parse schema JSON: {exc}", file=sys.stderr)
        return 1

    ok, msg = validate_with_jsonschema(instance, schema)
    if not ok:
        print(f"[validate][FAIL] {msg}", file=sys.stderr)
        return 1

    print(f"[validate][OK] {input_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
