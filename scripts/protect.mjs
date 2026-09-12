#!/usr/bin/env node
// Encrypts a markdown file so it can be embedded behind a client-side
// password gate via the {{< protected >}} shortcode.
//
// The password is never written anywhere — only the salt, IV, and
// ciphertext are printed. Those are safe to commit: without the password,
// PBKDF2 + AES-256-GCM makes brute-forcing them impractical for a
// reasonably strong password.
//
// Usage:
//   node scripts/protect.mjs "my strong password" path/to/secret.md
//
// Paste the printed shortcode into a content file in place of the
// plaintext you want hidden, then delete the plaintext scratch file.

import { randomBytes, pbkdf2Sync, createCipheriv } from "node:crypto";
import { readFileSync } from "node:fs";

const ITERATIONS = 250000;

const [, , password, inputPath] = process.argv;
if (!password || !inputPath) {
  console.error(
    "Usage: node scripts/protect.mjs <password> <path-to-plaintext-markdown-file>"
  );
  process.exit(1);
}

const plaintext = readFileSync(inputPath, "utf8");

const salt = randomBytes(16);
const iv = randomBytes(12);
const key = pbkdf2Sync(password, salt, ITERATIONS, 32, "sha256");

const cipher = createCipheriv("aes-256-gcm", key, iv);
const ciphertext = Buffer.concat([
  cipher.update(plaintext, "utf8"),
  cipher.final(),
]);
const authTag = cipher.getAuthTag();
const payload = Buffer.concat([ciphertext, authTag]);

const shortcode =
  `{{< protected salt="${salt.toString("base64")}" ` +
  `iv="${iv.toString("base64")}" ` +
  `iterations="${ITERATIONS}" ` +
  `data="${payload.toString("base64")}" >}}`;

console.log(shortcode);
