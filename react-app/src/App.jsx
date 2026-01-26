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

  let justClosed = false;

  for (const cmd of path.commands) {
    if (cmd.type === SVGPathData.MOVE_TO) {
      if (justClosed) {
        cubics.push({
          c0: { x, y: -y },
          c1: { x: (2 * x + cmd.x) / 3, y: (-2 * y -cmd.y) / 3 },
          c2: { x: (x + 2 * cmd.x) / 3, y: (-y - 2 * cmd.y) / 3 },
          c3: { x: cmd.x, y: -cmd.y },
        })
      } else {
        x = cmd.x
        y = cmd.y
      }
    } else if (cmd.type === SVGPathData.CURVE_TO) {
      cubics.push({
        c0: { x, y: -y },
        c1: { x: cmd.x1, y: -cmd.y1 },
        c2: { x: cmd.x2, y: -cmd.y2 },
        c3: { x: cmd.x, y: -cmd.y },
      })
      x = cmd.x
      y = cmd.y
    } else if (cmd.type === SVGPathData.LINE_TO) {
      cubics.push({
        c0: { x, y: -y },
        c1: { x: (2 * x + cmd.x) / 3, y: (-2 * y -cmd.y) / 3 },
        c2: { x: (x + 2 * cmd.x) / 3, y: (-y - 2 * cmd.y) / 3 },
        c3: { x: cmd.x, y: -cmd.y },
      })
      x = cmd.x
      y = cmd.y
    } else if (cmd.type == SVGPathData.CLOSE_PATH) {
      justClosed = true;
    } else {
      console.log("Not handling the following command")
      console.log(cmd)
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
    const v0 = getVertex(c.c0)
    const v1 = getVertex(c.c3)

    bs.add_edge(
      v0,
      v1,
      c.c1.x,
      c.c1.y,
      c.c2.x,
      c.c2.y
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
  const [target, setTarget] = useState(1)
  const [maxEdges, setMaxEdges] = useState(1)

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

  let bbox = {xmin:0, xmax:0, ymin:0, ymax:0};

  if (bs != null) {
    try {
      bbox = bs.get_bbox();
    } catch (e) {
      console.log(e.stack);
    }
  }

  const padding = 5

  const viewBox = [
    bbox.xmin - padding,
    bbox.ymin - padding,
    (bbox.xmax - bbox.xmin) + 2 * padding,
    (bbox.ymax - bbox.ymin) + 2 * padding,
  ].join(' ')

  let edges = [];

  if (bs != null) {
    try {
      edges = bs.edges();
    } catch (e) {
      console.log(e.stack);
    }
  }

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

        setMaxEdges(bs.number_of_edges())
        setBs(bs)
        forceUpdate()
      }}
    />

    {bs && (
      <input type="button" value="Initialize" onClick={(e) => { 
        if (bs) {
          bs.initialize()
          forceUpdate()
        }}}
        ></input>
    )}

    {bs && (
      <input type="button" value="Download SVG" onClick={(e) => {
        // Source - https://stackoverflow.com/a
        // Posted by DaveTheScientist, modified by community. See post 'Timeline' for change history
        // Retrieved 2026-01-26, License - CC BY-SA 4.0
        var svgData = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" + document.getElementById("graphSvg").outerHTML;
        var svgBlob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"});
        var svgUrl = URL.createObjectURL(svgBlob);
        var downloadLink = document.createElement("a");
        downloadLink.href = svgUrl;
        downloadLink.download = "output.svg";
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);
      }}></input>
    )}

    <br/>

    {bs && (
      <input
        type="range"
        min={0}
        max={maxEdges}
        step={1}
        value={target}
        onChange={(e) =>  {
          setTarget(Number(e.target.value))
          if (bs) {
            try {
              bs.run_to_complexity(e.target.value)
            } catch (e) {
              console.log(e.stack);
            }
          }
          forceUpdate()
        }}
      />
    )}

      {bs && (
        <svg version="1.2" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" id="graphSvg" width="100%" height={600} viewBox={viewBox}>
          {edges.map((eId) => {
            let curve;
            try {
              curve = bs.get_edge_curve(eId)
            } catch (e) {
              console.log(e.stack);
            }

            const d = `
              M ${curve.c0.x} ${curve.c0.y}
              C ${curve.c1.x} ${curve.c1.y},
                ${curve.c2.x} ${curve.c2.y},
                ${curve.c3.x} ${curve.c3.y}
            `

            return (
              <path
               	vector-effect="non-scaling-stroke"
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
