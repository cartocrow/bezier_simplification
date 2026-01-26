import { useEffect, useState, useReducer } from 'react'
import createModule from './wasm/wasm_frontend.js'
import { SVGPathData } from 'svg-pathdata'

function parseSVG(svgText) {
  const parser = new DOMParser()
  return parser.parseFromString(svgText, 'image/svg+xml')
}

function extractPaths(svgDoc) {
  return Array.from(svgDoc.querySelectorAll('path'))
    .map(p => p.getAttribute('d'))
    .filter(Boolean)
}

function pathToCubics(d) {
  const path = new SVGPathData(d)
    .toAbs()

  let x = 0, y = 0
  const cubics = []

  for (const cmd of path.commands) {
    if (cmd.type === SVGPathData.MOVE_TO) {
      x = cmd.x
      y = cmd.y
    }

    if (cmd.type === SVGPathData.CURVE_TO) {
      cubics.push({
        p0: { x, y },
        c0: { x: cmd.x1, y: cmd.y1 },
        c1: { x: cmd.x2, y: cmd.y2 },
        p1: { x: cmd.x, y: cmd.y },
      })
      x = cmd.x
      y = cmd.y
    }
  }

  return cubics
}

function addCubicsToGraph(bs, cubics) {
  const vertexMap = new Map()

  const getVertex = (pt) => {
    const key = `${pt.x},${pt.y}`
    if (!vertexMap.has(key)) {
      vertexMap.set(key, bs.insert_vertex(pt.x, pt.y))
    }
    return vertexMap.get(key)
  }

  for (const c of cubics) {
    const v0 = getVertex(c.p0)
    const v1 = getVertex(c.p1)

    bs.add_edge(
      v0,
      v1,
      c.c0.x,
      c.c0.y,
      c.c1.x,
      c.c1.y
    )
  }
}

function SvgDropZone({ onSvgLoaded }) {
  const onDrop = async (e) => {
    e.preventDefault()

    const file = e.dataTransfer.files[0]
    if (!file || !file.name.endsWith('.svg')) return

    const text = await file.text()
    onSvgLoaded(text)
  }

  return (
    <div
      onDragOver={(e) => e.preventDefault()}
      onDrop={onDrop}
      style={{
        border: '2px dashed #888',
        padding: '40px',
        textAlign: 'center',
        marginBottom: '1rem',
      }}
    >
      Drop SVG here
    </div>
  )
}


function App() {
  const [bs, setBs] = useState(null)

  const [ignored, forceUpdate] = useReducer(x => x + 1, 0)

  useEffect(() => {
    let cancelled = false
    let instance = null

    const init = async () => {
      const Module = await createModule({
        print: console.log,
        printErr: console.error,
      })

      if (cancelled) return

      console.log('WASM runtime initialized')

      const s = new Module.BezierSimplification()
      instance = s
      setBs(s)
    }

    init()

    return () => {
      cancelled = true
      instance?.delete()
    }
  }, [])

  const bbox = bs == null ? {xmin:0, xmax:0, ymin:0, ymax:0} : bs.get_bbox()

  const padding = 5

  const viewBox = [
    bbox.xmin - padding,
    bbox.ymin - padding,
    (bbox.xmax - bbox.xmin) + 2 * padding,
    (bbox.ymax - bbox.ymin) + 2 * padding,
  ].join(' ')

  return (
    <>
      <h1>Bezier simplification</h1>

    <SvgDropZone
      onSvgLoaded={(svgText) => {
        if (!bs) return

        const doc = parseSVG(svgText)
        const paths = extractPaths(doc)

        for (const d of paths) {
          const cubics = pathToCubics(d)
          addCubicsToGraph(bs, cubics)
        }

        setBs(bs)
        forceUpdate()
      }}
    />


      {bs && (
        <svg width="100%" height={800} viewBox={viewBox}>
          {bs.edges().map((eId) => {
            const curve = bs.get_edge_curve(eId)

            const d = `
              M ${curve.c0.x} ${curve.c0.y}
              C ${curve.c1.x} ${curve.c1.y},
                ${curve.c2.x} ${curve.c2.y},
                ${curve.c3.x} ${curve.c3.y}
            `

            return (
              <path
                key={eId}
                d={d}
                fill="none"
                stroke="black"
                strokeWidth={3}
              />
            )
          })}
          {/* {bs.vertices().map((vId) => {
            const pt = bs.get_vertex_point(vId)

            return (
              <circle
                key={vId}
                cx={pt.x}
                cy={pt.y}
                r={6}
                fill="white"
                stroke="black"
                strokeWidth={3}
              />
            )
          })} */}

        </svg>
      )}
    </>
  )
}

export default App
