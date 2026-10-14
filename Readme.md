# Personal Website

A static site built with [Hugo](https://gohugo.io), for quant/ML methodology
write-ups and experiment notes. No vendored theme — the whole look is one
plain stylesheet (`static/css/site.css`) and a handful of short templates
under `layouts/`. Open either one directly and edit it; there's no framework
indirection to fight.

**Design**: a sidebar-navigated "book" layout (persistent left nav, serif
headings, equation callout boxes, light/dark mode via CSS variables).
The sidebar does double duty: on most pages it shows the site's section tree;
on a page that's part of a `series` (see below), it shows that series' chapter
list instead, with the current chapter highlighted — i.e. it becomes the book's
table of contents.

## Requirements

- Hugo **extended** ≥ 0.146.0. Check with `hugo version`.
  - No package manager access? Download a binary from the
    [Hugo releases page](https://github.com/gohugoio/hugo/releases) — grab
    the `hugo_extended_<version>_linux-amd64.tar.gz` (or `darwin`/`windows`
    build for your OS), extract, and put the `hugo` binary on your `PATH`
    (e.g. `~/.local/bin`).

## Local development

```sh
./start.sh
```

Visit `http://localhost:1313`. Hugo live-reloads on file changes. The script
adds `~/.local/bin` to `PATH` (in case Hugo was installed there without
`sudo`) and forwards any extra flags to `hugo server`, e.g. `./start.sh
--port 1314`.

## Content structure

Posts live under section folders, not a flat `posts/`. The top-level
branches are:

- `content/overview/` — mindset and a high-level tour of the methodology,
  with the pipeline diagram. The entry point for a new reader.
- `content/quant/` — the quant research pipeline, one subfolder per stage:
  `universe/`, `feature-evaluation/`, `model/`, `alpha-evaluation/`.
- `content/machine-learning/` — general ML notes, not specific to the quant
  pipeline.

Each subfolder has an `_index.md` (section title + description) and zero or
more post files directly inside it. A subfolder's landing page
(`/quant/model/`, etc.) auto-lists whatever posts live inside it — no manual
index-page maintenance needed. `hugo.toml`'s `params.mainSections` (currently
`["quant", "machine-learning"]`) controls which top-level branches feed the
home page, archive, and search index.

The sidebar's site-tree (Overview / Quant stages / Machine Learning / About)
is hardcoded directly in `layouts/partials/sidebar.html` as a plain list of
`{slug, label}` pairs — not driven by `hugo.toml` menu config. If you add a
new top-level section, add one line there.

## Adding a post

```sh
hugo new quant/model/my-new-post.md
# or: hugo new machine-learning/my-new-post.md
```

Edit the front matter at the top of the file:

```yaml
---
title: "My New Post"
description: "One sentence — shown under the title and in listings."
date: 2026-09-12
tags: ["tag1", "tag2"]
math: true   # only needed if the post contains LaTeX
draft: false
---
```

Draft posts (`draft: true`) are excluded from `hugo build`/`hugo --minify`
but visible under `hugo server -D`.

### Math

Set `math: true` in a post's front matter to load KaTeX (see the
`{{ if .Params.math }}` block in `layouts/partials/head.html`). Then use
standard LaTeX delimiters in the body: `$inline$`, `$$block$$`,
`\(inline\)`, or `\[block\]`.

### Diagrams

`{{< mermaid >}}...{{< /mermaid >}}` renders any [Mermaid](https://mermaid.js.org/)
diagram definition. `{{< pipeline >}}` renders the
site's pipeline flowchart, with an optional comma-separated list of stages to
highlight (valid slugs: `universe`, `feature-evaluation`, `model`,
`alpha-evaluation`). Both load `mermaid.js` only on pages
that actually use them — see `layouts/shortcodes/{mermaid,pipeline}.html`
and the `.HasShortcode` check in `head.html`.

### Password-locking part of a post

Any section of any post can be hidden behind a password. It's encrypted
(AES-256-GCM, PBKDF2-derived key, 250k iterations) at authoring time and
decrypted entirely in the visitor's browser via the native Web Crypto API —
no backend, no build-pipeline changes, works on any static host. Without the
password, only ciphertext ever ships to the browser.

1. Write the content you want to hide in a scratch markdown file (anywhere
   outside `content/`, e.g. `/tmp/secret.md`).
2. Encrypt it:
   ```sh
   node scripts/protect.mjs "your password" /tmp/secret.md
   ```
   This prints a `{{< protected ... >}}` shortcode with the ciphertext baked
   in. The password itself is never written anywhere — only the salt, IV,
   and ciphertext are printed, and those are safe to commit.
3. Paste the printed shortcode into your post, in place of (or after) the
   plaintext. Delete the scratch file.
4. Delete your shell history entry for step 2 if it contains the real
   password (e.g. `history -d <line>`), since it was typed as a plain CLI
   argument.

The gated content supports the same Markdown + KaTeX math as the rest of the
site (rendered client-side via `marked.js` after a successful unlock — see
`layouts/shortcodes/protected.html` and the decrypt logic in
`layouts/partials/head.html`). If your protected content contains math,
also set `math: true` in the post's front matter as usual.

**Threat model, honestly stated:** this deters casual visitors and keeps the
plaintext out of "view source" / git history. It does not stop someone who
downloads the page and brute-forces the password offline — there's no rate
limiting on a static file. Don't use it for anything where a leak would
actually hurt; for real access control, gate the page at the edge instead
(e.g. Cloudflare Access or a Cloudflare Pages Function checking a password
before the page is ever served).

### Long-form articles and Table of Contents

Long write-ups are organized as unified, comprehensive guides (e.g., `multi-sub-universe-directional-coupling-geometry.md`, `cross-sectional-predictive-geometry.md`, `diffusion-statarb-implementation-guide.md`).

`layouts/partials/sidebar.html` automatically detects long articles and renders a sticky, scroll-spy Table of Contents (`Contents`) with quick-jump links and a "&larr; back to section" button, while preserving site-wide navigation below it.


### Drafts / source material

`docs/` holds raw notes that haven't been turned into posts yet. Converting
one is just: copy it into the right `content/quant/<stage>/` or
`content/machine-learning/` folder, add Hugo front matter, keep the math
delimiters as-is, and split into a series (see above) if it's long enough to
warrant it.

## Site configuration

Edit `hugo.toml` for `baseURL` (set this to your real domain or `*.pages.dev`
URL before deploying), `title`, `params.author`, `params.description`, and
`params.mainSections`. Edit `static/css/site.css` directly for anything
visual — colors are CSS custom properties at the top of the file (`--bg`,
`--ink`, `--feat`, etc.), duplicated once for the dark-mode media query and
once for the `[data-theme="dark"]` override (manual toggle).

## Deployment
 
### GitHub Pages (default)

Pushes to the `develop` branch trigger `.github/workflows/hugo.yml`, which automatically builds and deploys the site to GitHub Pages (`https://evanchen0830.github.io/`).

### Cloudflare Pages (alternative)

1. Push this repo to GitHub/GitLab.
2. In the Cloudflare dashboard: **Workers & Pages → Create → Pages → Connect to Git**, pick the repo.
3. Build settings:
   - Framework preset: **Hugo**
   - Build command: `hugo --gc --minify`
   - Build output directory: `public`
   - Environment variable: `HUGO_VERSION` = `0.150.1`
4. Deploy. Cloudflare rebuilds automatically on every push.

