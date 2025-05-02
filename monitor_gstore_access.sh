#!/bin/bash

WATCH_DIR="/srv/gstore/projects"
LOG_FILE="/shared/automount_access.log"

echo "$(date): Starting monitor as $(whoami) (UID=$(id -u), GID=$(id -g))" >> "$LOG_FILE"

# Ensure inotify-tools is installed
if ! command -v inotifywait >/dev/null 2>&1; then
  echo "inotify-tools is not installed. Exiting." >> "$LOG_FILE"
  exit 1
fi

inotifywait -m -e access "$WATCH_DIR" --format '%w%f' | while read path; do
  echo "$(date '+%F %T') Access detected on $path" >> "$LOG_FILE"
  ls "$path" > /dev/null 2>&1
done

