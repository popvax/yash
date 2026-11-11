// Public reads are limited to the exact verified catalog. The bucket stays private.
import objects from './objects.mjs';
import index from './index.mjs';
import manifest from './manifest.mjs';

function range(header, size) {
  if (!header) return null;
  const m = /^bytes=(\d*)-(\d*)$/.exec(header);
  if (!m || (!m[1] && !m[2]) || !size) return false;
  const start = m[1] ? Number(m[1]) : Math.max(0, size - Number(m[2]));
  const end = m[1] && m[2] ? Math.min(Number(m[2]), size - 1) : size - 1;
  if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start > end || start >= size) return false;
  return { offset: start, length: end - start + 1 };
}

export default {
  async fetch(request, env) {
    // The default deployment is inert. Public reads require an explicit launch.
    if (env.PUBLIC_DOWNLOADS_ENABLED !== 'true')
      return new Response('Downloads are not enabled', { status: 404, headers: { 'Cache-Control': 'no-store' } });
    if (!['GET', 'HEAD'].includes(request.method))
      return new Response('Read-only downloads', { status: 405, headers: { Allow: 'GET, HEAD' } });
    let key;
    try { key = decodeURIComponent(new URL(request.url).pathname.slice(1)); }
    catch { return new Response('Invalid path', { status: 400 }); }
    if (key === 'INDEX.json') return new Response(request.method === 'HEAD' ? null : JSON.stringify(index), {
      headers: { 'Content-Type': 'application/json', 'Cache-Control': 'no-store', 'X-Content-Type-Options': 'nosniff' },
    });
    if (key === index.object_manifest_key) return new Response(request.method === 'HEAD' ? null : manifest, {
      headers: {
        'Content-Type': 'application/json', 'Cache-Control': 'no-store',
        'Content-Length': String(new TextEncoder().encode(manifest).length),
        'X-Content-SHA256': index.object_manifest_sha256, 'X-Content-Type-Options': 'nosniff',
      },
    });
    if (!Object.hasOwn(objects, key)) return new Response(
      'Not served by this project. Obtain third-party tools, models and references from their providers: '+index.external_resources_url,
      { status: 404, headers: { 'Cache-Control': 'no-store', 'Link': `<${index.external_resources_url}>; rel="describedby"` } },
    );
    const expected = objects[key];
    const selected = range(request.headers.get('Range'), expected.bytes);
    if (selected === false) return new Response(null, { status: 416, headers: { 'Content-Range': `bytes */${expected.bytes}` } });
    const object = request.method === 'HEAD'
      ? await env.DATA.head(key)
      : await env.DATA.get(key, selected ? { range: selected } : undefined);
    if (!object) return new Response('Catalogued object unavailable', { status: 404 });
    if (object.size !== expected.bytes) return new Response('Catalog size mismatch', { status: 502 });
    const headers = new Headers({
      'Accept-Ranges': 'bytes',
      'Content-Type': 'application/octet-stream',
      'Content-Disposition': `attachment; filename="${key.split('/').pop().replace(/["\r\n]/g, '')}"`,
      'Content-Length': String(selected ? selected.length : object.size),
      'X-Content-Type-Options': 'nosniff',
      'X-Content-SHA256': expected.sha256,
      'ETag': object.httpEtag,
      'Cache-Control': 'public, max-age=3600',
      'Link': `<${index.third_party_notices_url}>; rel="describedby", <${index.external_resources_url}>; rel="describedby"`,
    });
    if (selected) headers.set('Content-Range', `bytes ${selected.offset}-${selected.offset + selected.length - 1}/${object.size}`);
    return new Response(request.method === 'HEAD' ? null : object.body, { status: selected ? 206 : 200, headers });
  },
};
