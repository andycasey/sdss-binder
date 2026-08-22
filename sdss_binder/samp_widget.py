"""
samp_widget — a SAMP Web Profile client that runs in the browser and delivers
messages to Python.

Because the JavaScript executes in *your* browser, the hub it talks to is the
one on *your* machine (http://localhost:21012). The Python process can live
anywhere — a remote marimo server, a container, a VM. No SSH tunnel required.

    import marimo as mo
    from samp_widget import SampClient

    client = SampClient(subscriptions=["sdss.spectrum.view"])
    w = mo.ui.anywidget(client)
    w

    # in another cell — re-runs every time TOPCAT sends a message
    msg = w.value["message"]
    sdss_id = msg["params"]["sdss_id"] if msg else None

Requires samp.js (the patched copy shipped alongside this file) in the same
directory.
"""

from __future__ import annotations

import base64
import pathlib
from typing import Any

import anywidget
import traitlets

_HERE = pathlib.Path(__file__).parent
_SAMP_JS = _HERE / "samp.js"


_WIDGET_JS = r"""
/* ------------------------------------------------------------------ *
 * anywidget front end. `samp` above is Mark Taylor's sampjs library.  *
 * ------------------------------------------------------------------ */

const BUILD = "3";

const PALETTE = {
  live:    "#3ecf8e",
  waiting: "#e0a33e",
  idle:    "#7c8b99",
  alert:   "#e5534b",
};

const CSS = `
.sampw {
  --sampw-line: color-mix(in srgb, currentColor 14%, transparent);
  --sampw-dim:  color-mix(in srgb, currentColor 55%, transparent);
  border: 1px solid var(--sampw-line);
  border-radius: 6px;
  padding: 10px 12px;
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  font-size: 13px;
  line-height: 1.45;
}
.sampw-head { display:flex; align-items:center; gap:10px; }
.sampw-dot {
  width:8px; height:8px; border-radius:50%; flex:none;
  background: ${PALETTE.idle};
  box-shadow: 0 0 0 0 transparent;
}
.sampw-dot[data-state="live"]    { background:${PALETTE.live};    animation: sampw-pulse 2.4s ease-out infinite; }
.sampw-dot[data-state="waiting"] { background:${PALETTE.waiting}; }
.sampw-dot[data-state="alert"]   { background:${PALETTE.alert};   }
@keyframes sampw-pulse {
  0%   { box-shadow: 0 0 0 0 color-mix(in srgb, ${PALETTE.live} 55%, transparent); }
  70%  { box-shadow: 0 0 0 6px transparent; }
  100% { box-shadow: 0 0 0 0 transparent; }
}
@media (prefers-reduced-motion: reduce) { .sampw-dot { animation: none !important; } }
.sampw-status { flex:1 1 auto; min-width:0; }
.sampw-count {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 11px; color: var(--sampw-dim); font-variant-numeric: tabular-nums;
}
.sampw-btn {
  font: inherit; font-size:12px; font-weight:600; line-height:1.2;
  padding:6px 14px; cursor:pointer; white-space:nowrap; flex:none;
  border-radius:5px; border:1px solid transparent;
  background:${PALETTE.live}; color:#08150f;
  box-shadow: 0 1px 2px rgba(0,0,0,.16);
  transition: background .12s ease, border-color .12s ease,
              box-shadow .12s ease, transform .06s ease;
}
.sampw-btn:hover:not(:disabled) {
  background: color-mix(in srgb, ${PALETTE.live} 82%, white);
  box-shadow: 0 2px 5px rgba(0,0,0,.2);
}
.sampw-btn:active:not(:disabled) { transform: translateY(1px); box-shadow:none; }
.sampw-btn:focus-visible { outline:2px solid ${PALETTE.live}; outline-offset:2px; }
/* Once connected, disconnecting is the rare action - step it back to an outline. */
.sampw-btn[data-variant="secondary"] {
  background:transparent; color:inherit; box-shadow:none; font-weight:500;
  border-color: color-mix(in srgb, currentColor 32%, transparent);
}
.sampw-btn[data-variant="secondary"]:hover:not(:disabled) {
  background: color-mix(in srgb, currentColor 8%, transparent);
  border-color: color-mix(in srgb, currentColor 55%, transparent);
  box-shadow:none;
}
.sampw-btn:disabled {
  background:transparent; color: var(--sampw-dim); font-weight:500;
  border-color: var(--sampw-line); box-shadow:none; cursor:not-allowed;
}
@media (prefers-reduced-motion: reduce) { .sampw-btn { transition:none; } }
.sampw-log {
  margin-top:9px; padding-top:8px; border-top:1px solid var(--sampw-line);
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size:11.5px;
}
.sampw-log:empty { display:none; }
.sampw-row { display:flex; gap:8px; padding:1px 0; }
.sampw-t { color: var(--sampw-dim); font-variant-numeric: tabular-nums; flex:none; }
.sampw-m { flex:none; }
.sampw-p { color: var(--sampw-dim); overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
.sampw-hint { color: var(--sampw-dim); font-size:12px; margin-top:6px; }
.sampw-hint:empty { display:none; }
`;

function b64(buf) {
  const bytes = new Uint8Array(buf);
  let out = "";
  for (let i = 0; i < bytes.length; i += 0x8000) {
    out += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
  }
  return btoa(out);
}

function brief(params) {
  return Object.entries(params || {})
    .map(([k, v]) => k + "=" + String(v).slice(0, 48))
    .join("  ");
}

function render({ model, el }) {
  // The stylesheet has to live in the same tree as `el`. Some hosts render
  // widgets inside a shadow root, where a <style> appended to document.head
  // never reaches us — only inherited properties like font and color cross the
  // boundary, so the widget renders as unstyled text. Appending inside `el`
  // itself is correct in both light DOM and shadow DOM. It is also rewritten on
  // every render, since marimo re-runs cells without reloading the page.
  let style = el.querySelector("style#sampw-css");
  if (!style) {
    style = document.createElement("style");
    style.id = "sampw-css";
    el.appendChild(style);
  }
  style.textContent = CSS;
  style.dataset.build = BUILD;

  // Earlier builds appended their stylesheet to document.head. Clear any that
  // are still around so they cannot compete with the rules above.
  document.head
    .querySelectorAll("style#sampw-css")
    .forEach((stale) => stale.remove());

  const root = document.createElement("div");
  root.className = "sampw";
  root.dataset.build = BUILD;
  root.innerHTML = `
    <div class="sampw-head">
      <span class="sampw-dot" data-state="idle"></span>
      <span class="sampw-status"></span>
      <span class="sampw-count"></span>
      <button class="sampw-btn" type="button"></button>
    </div>
    <div class="sampw-hint"></div>
    <div class="sampw-log"></div>`;
  el.appendChild(root);

  const $dot    = root.querySelector(".sampw-dot");
  const $status = root.querySelector(".sampw-status");
  const $count  = root.querySelector(".sampw-count");
  const $btn    = root.querySelector(".sampw-btn");
  const $hint   = root.querySelector(".sampw-hint");
  const $log    = root.querySelector(".sampw-log");

  let connection = null;
  let hubUp = false;
  let seq = model.get("message_count") || 0;
  const history = [];

  const paint = () => {
    const registered = !!connection;
    $dot.dataset.state = registered ? "live" : hubUp ? "waiting" : "idle";
    if (registered) {
      $status.textContent = "Listening as " + model.get("client_name");
    } else if (hubUp) {
      $status.textContent = "Hub found on this machine";
    } else {
      $status.textContent = "No hub on this machine";
    }
    $count.textContent = seq ? seq + (seq === 1 ? " message" : " messages") : "";
    $btn.textContent = registered ? "Disconnect" : "Connect to hub";
    $btn.dataset.variant = registered ? "secondary" : "primary";
    $btn.disabled = !hubUp && !registered;
    $btn.title = registered
      ? "Unregister this page from the SAMP hub"
      : hubUp
      ? "Register this page with the SAMP hub running on your machine"
      : "No SAMP hub is running on your machine";
    $hint.textContent = registered
      ? ""
      : hubUp
      ? "TOPCAT will ask you to approve this page before it can receive anything."
      : "Start TOPCAT (or another SAMP hub) on the computer running this browser.";
  };

  const drawLog = () => {
    $log.innerHTML = history
      .slice(-5)
      .reverse()
      .map(
        (h) =>
          `<div class="sampw-row"><span class="sampw-t">${h.t}</span>` +
          `<span class="sampw-m">${h.mtype}</span>` +
          `<span class="sampw-p">${h.brief}</span></div>`
      )
      .join("");
  };

  /* ---------------- receive: hub -> browser -> python ---------------- */

  const tracker = new samp.ClientTracker();
  const maxBytes = model.get("max_fetch_bytes");

  const deliver = async (senderId, message, isCall) => {
    const mtype = message["samp.mtype"];
    const params = message["samp.params"] || {};
    const meta = tracker.metas[senderId] || {};
    const payload = {
      seq: ++seq,
      mtype: mtype,
      sender: senderId,
      sender_name: meta["samp.name"] || senderId,
      params: params,
      received: new Date().toISOString(),
    };

    // A file:// URL from TOPCAT is meaningless to a remote Python process and
    // unfetchable from a web page. The hub exposes a translator that serves it
    // over HTTP to this origin; hand Python both forms.
    if (params.url && connection) {
      try {
        payload.url_translated = connection.translateUrl(params.url);
      } catch (e) {
        /* translator absent; leave unset */
      }
    }

    if (maxBytes > 0 && payload.url_translated && /^(table|image)\.load\./.test(mtype)) {
      try {
        const r = await fetch(payload.url_translated);
        if (!r.ok) throw new Error("HTTP " + r.status);
        const buf = await r.arrayBuffer();
        if (buf.byteLength > maxBytes) {
          payload.fetch_error =
            "file is " + buf.byteLength + " bytes, over the " + maxBytes + " byte limit";
        } else {
          payload.data_b64 = b64(buf);
          payload.data_bytes = buf.byteLength;
        }
      } catch (e) {
        payload.fetch_error = String(e);
      }
    }

    model.set("message", payload);
    model.set("message_count", payload.seq);
    model.save_changes();

    history.push({
      t: new Date().toLocaleTimeString([], { hour12: false }),
      mtype: mtype,
      brief: brief(params),
    });
    drawLog();
    paint();
  };

  const wanted = (model.get("subscriptions") || []).filter(
    (mt) => !mt.startsWith("samp.hub.")
  );
  for (const mt of wanted) tracker.callHandler[mt] = deliver;

  // sampjs dispatches on an exact mtype key, so a wildcard subscription would
  // be declared to the hub but never routed. Bind the concrete mtype the first
  // time we see one that matches a wildcard we asked for.
  const subs = tracker.calculateSubscriptions();
  const baseNotify = tracker.receiveNotification.bind(tracker);
  const baseCall = tracker.receiveCall.bind(tracker);
  const bindWildcard = (mtype) => {
    if (mtype in tracker.callHandler) return;
    for (const pat of wanted) {
      if (pat.indexOf("*") >= 0 && samp.isSubscribed(subs, mtype)) {
        tracker.callHandler[mtype] = deliver;
        return;
      }
    }
  };
  tracker.receiveNotification = (senderId, message) => {
    bindWildcard(message["samp.mtype"]);
    return baseNotify(senderId, message);
  };
  tracker.receiveCall = (senderId, msgId, message) => {
    bindWildcard(message["samp.mtype"]);
    return baseCall(senderId, msgId, message);
  };

  const meta = {
    "samp.name": model.get("client_name"),
    "samp.description.text": model.get("client_description"),
  };
  if (model.get("icon_url")) meta["samp.icon.url"] = model.get("icon_url");

  const connector = new samp.Connector(model.get("client_name"), meta, tracker, subs);

  connector.onreg = (conn) => {
    connection = conn;
    model.set("registered", true);
    model.set("status", "registered");
    model.save_changes();
    paint();
  };
  connector.onunreg = () => {
    connection = null;
    model.set("registered", false);
    model.set("status", "unregistered");
    model.save_changes();
    paint();
  };

  $btn.onclick = () => {
    if (connection) connector.unregister();
    else connector.register();
  };

  const timer = connector.onHubAvailability((up) => {
    if (up === hubUp) return;
    hubUp = up;
    model.set("hub_available", up);
    model.save_changes();
    if (up && model.get("autoconnect") && !connection) connector.register();
    paint();
  }, 2000);

  /* ---------------- send: python -> browser -> hub ---------------- */

  const onOutbox = () => {
    const out = model.get("outbox");
    if (!out || !out.mtype) return;
    const send = (conn) => {
      const msg = new samp.Message(out.mtype, out.params || {});
      if (out.target) conn.notify([out.target, msg]);
      else conn.notifyAll([msg]);
    };
    if (connection) send(connection);
    else connector.runWithConnection(send);
  };
  model.on("change:outbox", onOutbox);

  paint();

  return () => {
    clearInterval(timer);
    model.off("change:outbox", onOutbox);
    try {
      connector.unregister();
    } catch (e) {
      /* page is going away anyway */
    }
  };
}

export default { render };
"""


def _build_esm() -> str:
    if not _SAMP_JS.exists():
        raise FileNotFoundError(
            f"samp.js not found at {_SAMP_JS}. Keep it next to samp_widget.py."
        )
    lib = _SAMP_JS.read_text(encoding="utf-8")
    missing = [
        name
        for marker, name in (
            ("var TYPE_STRING", "TYPE_STRING/TYPE_LIST/TYPE_MAP"),
            ("var oc = this.onclose", "oc in Connection.close"),
        )
        if marker not in lib
    ]
    if missing:
        raise ValueError(
            "This samp.js still assigns "
            + " and ".join(missing)
            + " as implicit globals, which throws under the strict mode that ES "
            "modules always run in. Use the patched copy shipped alongside this file."
        )
    return lib + "\n" + _WIDGET_JS


class SampClient(anywidget.AnyWidget):
    """A SAMP Web Profile client living in the browser, reporting to Python."""

    _esm = _build_esm()

    # configuration
    client_name = traitlets.Unicode("marimo").tag(sync=True)
    client_description = traitlets.Unicode("Notebook SAMP client").tag(sync=True)
    icon_url = traitlets.Unicode("").tag(sync=True)
    subscriptions = traitlets.List(traitlets.Unicode(), default_value=[]).tag(sync=True)
    autoconnect = traitlets.Bool(False).tag(sync=True)
    # Fetched table bytes travel over the notebook websocket as base64, so keep
    # this modest. 0 disables fetching and passes only URLs through.
    max_fetch_bytes = traitlets.Int(0).tag(sync=True)

    # live state
    hub_available = traitlets.Bool(False).tag(sync=True)
    registered = traitlets.Bool(False).tag(sync=True)
    status = traitlets.Unicode("").tag(sync=True)

    # inbound
    message = traitlets.Dict(allow_none=True, default_value=None).tag(sync=True)
    message_count = traitlets.Int(0).tag(sync=True)

    # outbound
    outbox = traitlets.Dict(default_value={}).tag(sync=True)

    def __init__(self, subscriptions: list[str] | None = None, **kwargs: Any) -> None:
        if subscriptions is not None:
            kwargs["subscriptions"] = list(subscriptions)
        super().__init__(**kwargs)
        self._outbox_seq = 0

    def send(
        self,
        mtype: str,
        params: dict[str, Any] | None = None,
        target: str = "",
    ) -> None:
        """Send a SAMP message out through the hub.

        target is a client id; the default broadcasts to every subscriber.
        Values must be SAMP-legal: strings, or lists/dicts thereof.
        """
        self._outbox_seq += 1
        self.outbox = {
            "seq": self._outbox_seq,
            "mtype": mtype,
            "params": params or {},
            "target": target,
        }

    def data(self) -> bytes | None:
        """Bytes of the table attached to the last message, if any were fetched."""
        msg = self.message or {}
        blob = msg.get("data_b64")
        return base64.b64decode(blob) if blob else None
