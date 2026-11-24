// Cloudflare Pages: parola d'ordine unica (HTTP Basic) davanti a tutto il sito.
// Variabile d'ambiente da impostare in Pages → Settings → Variables and Secrets:
//   PIGATO_PASSWORD  (secret)
// Nome utente: qualsiasi, oppure PIGATO_UTENTE se impostato.

function uguali(a, b) {
  // confronto a tempo costante
  const x = new TextEncoder().encode(a), y = new TextEncoder().encode(b);
  let diff = x.length ^ y.length;
  for (let i = 0; i < Math.max(x.length, y.length); i++) diff |= (x[i] ?? 0) ^ (y[i] ?? 0);
  return diff === 0;
}

function negato(msg = "Accesso riservato") {
  return new Response(msg, {
    status: 401,
    headers: { "WWW-Authenticate": 'Basic realm="PIGATO TDM", charset="UTF-8"', "Cache-Control": "no-store" },
  });
}

export async function onRequest({ request, env, next }) {
  const attesa = env.PIGATO_PASSWORD;
  if (!attesa) return new Response("Parola d'ordine non configurata (PIGATO_PASSWORD).", { status: 503 });
  const h = request.headers.get("Authorization") || "";
  if (!h.startsWith("Basic ")) return negato();
  let utente = "", pw = "";
  try {
    const dec = new TextDecoder().decode(Uint8Array.from(atob(h.slice(6)), c => c.charCodeAt(0)));
    const i = dec.indexOf(":");
    utente = dec.slice(0, i); pw = dec.slice(i + 1);
  } catch { return negato(); }
  if (!uguali(pw, attesa) || (env.PIGATO_UTENTE && !uguali(utente, env.PIGATO_UTENTE))) return negato();
  const risposta = await next();
  const r = new Response(risposta.body, risposta);
  r.headers.set("X-Robots-Tag", "noindex, nofollow");
  r.headers.set("Cache-Control", "private, no-store");
  return r;
}
