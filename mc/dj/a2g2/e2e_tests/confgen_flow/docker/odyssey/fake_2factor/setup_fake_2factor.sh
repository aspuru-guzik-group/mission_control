#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Compile pam module for fake 2factor auth.
gcc \
  -fPIC \
  -fno-stack-protector \
  -c $DIR/fake_2factor_pam_module.c \
  -o /fake_2factor.o

ld \
  -x \
  --shared \
  -o /lib64/security/pam_fake_2factor.so \
  /fake_2factor.o

# Configure pam auth for sshd to use fake 2factor module.
FAKE_2FACTOR_PAM_AUTH_LINE="auth required pam_fake_2factor.so"
sed -i "/auth.*substack.*password-auth/a $FAKE_2FACTOR_PAM_AUTH_LINE" /etc/pam.d/sshd

# Configure sshd to use 2factor, via challengeresponse.
sed -i "s/ChallengeResponseAuthentication no/ChallengeResponseAuthentication yes/" /etc/ssh/sshd_config
