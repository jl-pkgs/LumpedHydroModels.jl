#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.2": plot

#import "cetz.typ": *

#show "{": ""
#show "}": ""

#let b = 0.4
#let cal_MM(x) = 1 - calc.pow(1 - x, 1 / b)
#let cal_a(MM) = 1 - calc.pow(1 - MM, b)

#let a0 = 0.2
#let MM0 = cal_MM(a0)
#let MM1 = 0.65
#let a1 = cal_a(MM1)

#let MM = "MM'"
#let WM = "WM'"
#let WMM = "WMM"
#let dWM = $d"WM'"$

#let ifelse(con, a, b) = if (con) { a } else { b }

#let h_arrow(width, x0, y0, name: "Q0", text: none, arrow: "straight", anchor: "west", ..args) = {
  line(
    (x0, y0),
    (x0 + width, y0),
    mark: (end: "stealth", fill: black),
    name: name,
    ..args,
  )
  content(name + ".end", label(outset: -1pt, inset: 3pt)[#text], anchor: anchor)
}

#let v_arrow(height, x0, y0, name: "Q0", text: none, arrow: "straight", anchor: "west", ..args) = {
  line(
    (x0, y0),
    (x0, y0 - height),
    mark: (end: "stealth", fill: black),
    name: name,
    ..args,
  )
  content(name + ".end", label(outset: -1pt, inset: 3pt)[#text], anchor: anchor)
}

#let axis-config = (
  size: (1, 1),
  axis-style: "left",
  x-label: [$bold(alpha)$], // Y轴标题
  x-tick-step: none,
  y-tick-step: none,
  y-min: 0,
  y-max: 1,
  x-min: 0,
  x-max: 1,
)

#let label-title(title) = {
  content((-0.15, 1.01), text(fill: black, size: 13pt)[#title])
}


#let fig-W = cetz.canvas(length: 50mm, {
  set-style(stroke: (thickness: 0.8pt))
  set-style(
    axes: (stroke: 0.7pt, tick: (stroke: .5pt)),
    legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 80%),
  )

  plot.plot(
    y-label: [$bold("WM")$], // X轴标题
    ..axis-config,
    {
      plot.add(cal_MM, domain: (0, 1))
      // W
      plot.add-fill-between(
        domain: (0, 1),
        style: (fill: black.lighten(90%), stroke: none),
        x => ifelse(x > a0, MM0, cal_MM(x)),
        x => { 0 },
      )
      // dW
      plot.add-fill-between(
        domain: (a0, 1),
        style: (fill: blue.lighten(80%), stroke: none),
        x => ifelse(x > a1, MM1, cal_MM(x)),
        x => { MM0 },
      )
      // dR
      plot.add-fill-between(
        domain: (0, a1),
        style: (fill: green.lighten(80%), stroke: none),
        x => { MM1 },
        x => ifelse(x < a0, MM0, cal_MM(x)),
      )
    },
  )
  // y_axis(1.1, 0, 0, name: "y", text: "WM'", anchor: "south")
  // x_axis(1.1, 0, 0, name: "x", text: $alpha$)

  hline(1, 1, stroke: (dash: "dashed", paint: black.lighten(60%)))
  vline(1, 1, stroke: (dash: "dashed", paint: black.lighten(60%)))

  line((0, MM0), (1, MM0), stroke: (dash: "dashed", paint: blue.lighten(40%)))
  line((0, MM1), (1, MM1), stroke: (dash: "solid", paint: blue.lighten(40%)))
  content((0.6, 0.2), text(fill: black)[$bold(W)$])

  content((0.6, (MM0 + MM1) / 2), text(fill: blue)[$bold(Delta W)$])

  content((0.5, 1.07), $bold(alpha) = 1 - (1 - WM/WMM)^b$)

  content((a1 / 2 - 0.05, (MM0 + MM1) / 2), text(fill: green)[$bold("R")$])

  mv_arrow(MM1 - MM0, -0.04, MM1, text: $"PE"$)
  mv_arrow(MM0 - 0, -0.04, MM0, text: $a$)
  // content((1, -0.03), label($1.0$), anchor: "north")
  label-title($bold((a))$) // title
})

#let fig-S = cetz.canvas(length: 50mm, {
  set-style(stroke: (thickness: 0.8pt))
  set-style(
    axes: (stroke: 0.7pt, tick: (stroke: .5pt)),
    legend: (stroke: none, orientation: ttb, item: (spacing: .3), scale: 80%),
  )
  
  plot.plot(
    y-label: [$bold("S")$], // X轴标题
    ..axis-config,
    {
      plot.add(cal_MM, domain: (0, 1))
      // W
      plot.add-fill-between(
        domain: (0, 1),
        style: (fill: black.lighten(90%), stroke: none),
        x => ifelse(x > a0, MM0, cal_MM(x)),
        x => { 0 },
      )
      // dW
      plot.add-fill-between(
        domain: (a0, 1),
        style: (fill: blue.lighten(80%), stroke: none),
        x => ifelse(x > a1, MM1, cal_MM(x)),
        x => { MM0 },
      )
      // dR
      plot.add-fill-between(
        domain: (0, a1),
        style: (fill: green.lighten(80%), stroke: none),
        x => { MM1 },
        x => ifelse(x < a0, MM0, cal_MM(x)),
      )
    },
  )
  
  hline(1, 1, stroke: (dash: "dashed", paint: black.lighten(60%)))
  vline(1, 1, stroke: (dash: "dashed", paint: black.lighten(60%)))
  
  line((0, MM0), (1, MM0), stroke: (dash: "dashed", paint: blue.lighten(40%)))
  line((0, MM1), (1, MM1), stroke: (dash: "solid", paint: blue.lighten(40%)))
  
  content((0.6, 0.2), text(fill: black)[$bold(S_1)$])
  content((0.6, (MM0 + MM1) / 2), text(fill: blue)[$bold(Delta S)$])
  
  content((0.5, 1.07), $bold(alpha) = 1 - (1 - S/S_"mm")^"EX"$)
  content((a1 / 2 - 0.05, (MM0 + MM1) / 2), text(fill: green)[$bold("RS")$])
  
  let mar = 0.035
  v_arrow(0.1, 0.6, mar, name: "RG", text: $bold("RG")$, anchor: "north")
  h_arrow(0.1, 1.0 - mar, name: "RI", text: $bold("RI")$, 0.1)
  
  mv_arrow(MM1 - MM0, -0.04, MM1, text: $"PE"$)
  mv_arrow(MM0 - 0, -0.04, MM0, text: $"AU"$)

  label-title($bold((b))$) // title
})

#figure(
  grid(
    columns: 2,
  )[#fig-W][#fig-S],
  caption: [
    新安江模型蓄水容量曲线（a）与自由水蓄量曲线（b）。
  ],
) <fig_>
