# Personal Website

A static blog built with [Hugo](https://gohugo.io) and the
[PaperMod](https://github.com/adityatelange/hugo-PaperMod) theme, for
methodology write-ups and experiment notes. The theme is vendored directly
into `themes/PaperMod` (plain files, no git submodule), so `git clone` alone
gets you a buildable site — nothing extra to init.

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

Posts live under section folders, not a flat `posts/`. The two top-level
branches are:

- `content/quant/` — the quant research pipeline, one subfolder per stage,
  in the order a signal actually moves through it:
  `universe/`, `features/`, `feature-selection/`, `model/`, `post-alpha/`,
  `portfolio/`.
- `content/machine-learning/` — general ML notes, not specific to the quant
  pipeline.

Each subfolder has an `_index.md` (section title + description) and zero or
more post files directly inside it. A subfolder's landing page
(`/quant/model/`, etc.) auto-lists whatever posts live inside it — no manual
index-page maintenance needed. `hugo.toml`'s `params.mainSections` (currently
`["quant", "machine-learning"]`) controls which top-level branches feed the
home page, archive, and RSS feed.

## Adding a post

```sh
hugo new quant/model/my-new-post.md
# or: hugo new machine-learning/my-new-post.md
```

Edit the front matter at the top of the file:

```yaml
---
title: "My New Post"
date: 2026-09-12
tags: ["tag1", "tag2"]
math: true   # only needed if the post contains LaTeX
draft: false
---
```

Draft posts (`draft: true`) are excluded from `hugo build`/`hugo --minify`
but visible under `hugo server -D`.

### Math

Set `math: true` in a post's front matter to load KaTeX (see
`layouts/partials/extend_head.html`). Then use standard LaTeX delimiters in
the body: `$inline$`, `$$block$$`, `\(inline\)`, or `\[block\]`.

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
`layouts/partials/extend_head.html`). If your protected content contains
math, also set `math: true` in the post's front matter as usual.

**Threat model, honestly stated:** this deters casual visitors and keeps the
plaintext out of "view source" / git history. It does not stop someone who
downloads the page and brute-forces the password offline — there's no rate
limiting on a static file. Don't use it for anything where a leak would
actually hurt; for real access control, gate the page at the edge instead
(e.g. Cloudflare Access or a Cloudflare Pages Function checking a password
before the page is ever served).

### Multi-part series

To turn a long write-up into a linked series (like
`content/quant/feature-selection/part-1-*.md` through `part-4-*.md`), give
each post a shared `series` name and a `part` number in front matter:

```yaml
series: ["My Series Name"]
part: 1
```

Then drop `{{< series name="My Series Name" >}}` right after the front
matter in each post's body. It renders a small box linking every post that
shares that series name, ordered by `part`, bolding whichever one you're
currently reading. The shortcode lives at
`layouts/shortcodes/series.html` — it's plain Hugo, not theme-specific, so it
survives theme updates.

### Drafts / source material

`docs/` holds raw notes that haven't been turned into posts yet. Converting
one is just: copy it into the right `content/quant/<stage>/` or
`content/machine-learning/` folder, add Hugo front matter, keep the math
delimiters as-is, and split into a series (see above) if it's long enough to
warrant it.

## Site configuration

Edit `hugo.toml`:

- `baseURL` — set this to your real domain (or `*.pages.dev` URL) before
  deploying, so canonical links, RSS, and the sitemap are correct.
- `title`, `params.author`, `params.description` — site identity.
- `params.homeInfoParams` — the intro blurb shown on the home page.
- `params.socialIcons` — links shown in the header/footer (GitHub, email,
  RSS, etc. — see the
  [PaperMod wiki](https://github.com/adityatelange/hugo-PaperMod/wiki) for
  the full list of supported icon names).
- `menu.main` — top nav items.

## Deploying to Cloudflare Pages

1. Push this repo to GitHub/GitLab.
2. In the Cloudflare dashboard: **Workers & Pages → Create → Pages → Connect
   to Git**, pick the repo.
3. Build settings:
   - Framework preset: **Hugo**
   - Build command: `hugo --gc --minify`
   - Build output directory: `public`
   - Environment variable: `HUGO_VERSION` = `0.150.1` (or whatever
     `hugo version` prints locally — keep them in sync)
4. Deploy. Cloudflare rebuilds automatically on every push.
5. To use a custom domain, add it under the Pages project's **Custom
   domains** tab; otherwise the site is live at `<project>.pages.dev`.

## Updating the theme

The theme is a plain vendored copy, not a submodule, so it won't drift or
require `git submodule update`. To pull in upstream PaperMod changes later:

```sh
git clone --depth 1 https://github.com/adityatelange/hugo-PaperMod.git /tmp/papermod
rsync -a --exclude='.git' /tmp/papermod/ themes/PaperMod/
```

Then re-check `layouts/partials/extend_head.html` and any other files under
the site's own `layouts/` still apply cleanly (they override the theme's
files of the same name, so upstream changes to those specific files won't
take effect until you diff and merge manually).
