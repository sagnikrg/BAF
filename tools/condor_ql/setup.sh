#!/bin/bash
# setup.sh: installs condor_ql + batch_time_left for the current BAF user (safe to run again; it also updates).
#   curl -fsSL https://raw.githubusercontent.com/sagnikrg/BAF/main/tools/condor_ql/setup.sh | bash
set -e
DIR="$HOME/.local/condor_ql"
URL="https://raw.githubusercontent.com/sagnikrg/BAF/main/tools/condor_ql"
mkdir -p "$DIR"
for f in condor_ql batch_time_left LICENSE README.md; do
  curl -fsSL "$URL/$f" -o "$DIR/$f"
done
chmod +x "$DIR/condor_ql" "$DIR/batch_time_left"

# PATH + "cd ~ first" wrapper in ~/.bashrc, once
if ! grep -q "# condor_ql setup" ~/.bashrc 2>/dev/null; then
  cat >> ~/.bashrc <<'EOF'

# condor_ql setup (BAF/tools/condor_ql)
export PATH="$HOME/.local/condor_ql:$PATH"
condor_ql() { cd ~ && command condor_ql "$@"; }
EOF
fi

# ssh logins read ~/.bash_profile (if it exists) instead of ~/.profile: make sure it loads ~/.bashrc
if [ -f ~/.bash_profile ] && ! grep -q "bashrc" ~/.bash_profile; then
  echo '[ -f ~/.bashrc ] && . ~/.bashrc' >> ~/.bash_profile
fi

echo "condor_ql installed in $DIR. Open a new ssh session (or: source ~/.bashrc), then run: condor_ql"
