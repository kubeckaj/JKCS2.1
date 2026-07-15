import os
import sys
import threading
from datetime import datetime

# Process-wide locks so concurrent monitoring threads cannot tear lines.
_TERM_LOCK = threading.Lock()
_FILE_LOCK = threading.Lock()

# level -> (log file tag, terminal prefix)
_LEVELS = {
    'info':    ('INFO ', ''),
    'event':   ('INFO ', ''),
    'success': ('OK   ', ''),
    'warning': ('WARN ', 'WARNING: '),
    'error':   ('ERROR', 'ERROR: '),
}


def _timestamp():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def _write_terminal(stream, text):
    with _TERM_LOCK:
        stream.write(text + '\n')
        stream.flush()


class Logger:
    def __init__(self, log_file=None, tag=None):
        self.bind(log_file, tag)

    def bind(self, log_file, tag=None):
        self.log_file = log_file
        if tag is None and log_file:
            tag = os.path.basename(os.path.dirname(os.path.abspath(log_file)))
        self.tag = tag

    def _write_file(self, level_tag, message):
        if not self.log_file:
            return
        stamp = _timestamp()
        lines = ''.join(f"{stamp} [{level_tag}] {line}\n" for line in message.split('\n'))
        with _FILE_LOCK:
            with open(self.log_file, 'a') as f:
                f.write(lines)

    def _emit(self, level, message, echo):
        file_tag, term_prefix = _LEVELS[level]
        self._write_file(file_tag, message)
        if not echo:
            return
        stream = sys.stderr if level == 'error' else sys.stdout
        tag = f"[{self.tag}] " if self.tag else ''
        _write_terminal(stream, f"{tag}{term_prefix}{message}")

    def info(self, message):
        self._emit('info', message, echo=self.log_file is None)

    def event(self, message):
        self._emit('event', message, echo=True)

    def success(self, message):
        self._emit('success', message, echo=True)

    def warning(self, message):
        self._emit('warning', message, echo=True)

    def error(self, message):
        self._emit('error', message, echo=True)

    def results(self, text):
        self._write_file('OK   ', text)
        _write_terminal(sys.stdout, box(text))


# The one logger of this process. main() binds it to the working directory's log
# file; until then (and for runs that own no directory) it writes to the terminal.
logger = Logger()


def box(text):
    lines = text.split('\n')
    rule = '=' * (max(len(line) for line in lines) + 2)
    body = '\n'.join(f" {line}" for line in lines)
    return f"{rule}\n{body}\n{rule}"


def banner(pairs, title='JKTS'):
    label_width = max(len(label) for label, _ in pairs)
    lines = [f" {label.ljust(label_width)}   {value}" for label, value in pairs]
    rule = '=' * max(len(line) for line in lines)
    _write_terminal(sys.stdout, '\n'.join([title, rule] + lines + [rule]))
