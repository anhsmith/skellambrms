# Generate the skellambrms hex sticker: man/figures/logo.svg and logo.png.
# Run from the package root.
#
# The motif is the model itself. What is left after two paired counts are
# differenced is a single Z-valued distribution centred near zero: a low-count
# Skellam pmf, drawn as a bar for each integer in -4:4. Grey marks exact
# agreement (d = 0), teal the mass on the positive side, orange the mass on the
# negative side. The palette matches the companion package bicountbrms, where
# the same two colours mark the two sources' private components.

W <- 519.6   # hex width  (1.732 : 2, the hexb.in ratio)
H <- 600     # hex height

col_zero <- "#B4AA98"  # warm grey -- d = 0, the two sources agree
col_pos  <- "#5BA5C7"  # teal      -- d > 0, source 1 counts higher
col_neg  <- "#E28A45"  # orange    -- d < 0, source 2 counts higher
col_ink  <- "#2F5F7A"
col_bg   <- "#FCFAF6"

hex_pts <- paste(
  sprintf(
    "%.2f,%.2f",
    c(W / 2, W, W, W / 2, 0, 0),
    c(0, H / 4, 3 * H / 4, H, 3 * H / 4, H / 4)
  ),
  collapse = " "
)

# ---- the pmf ----------------------------------------------------------------
# a symmetric, low-count Skellam: mu1 = mu2 = 2, so four bars either side of
# zero before the mass is negligible

d   <- -4:4
pmf <- skellam::dskellam(d, lambda1 = 2, lambda2 = 2)

bar_w   <- 32
gap_x   <- 8
pitch_x <- bar_w + gap_x
h_max   <- 228                       # height of the tallest bar
h_min   <- 13                        # floor, so the outermost bars stay visible

heights <- pmax(h_max * pmf / max(pmf), h_min)

base_y <- 340                                        # common baseline
span   <- length(d) * pitch_x - gap_x
x0     <- (W - span) / 2

# bars are rounded at the top only, so they sit flat on the axis
bar_r <- 5

bars <- character(0)

for (i in seq_along(d)) {
  fill <- if (d[i] == 0) col_zero else if (d[i] > 0) col_pos else col_neg
  x    <- x0 + (i - 1) * pitch_x
  y    <- base_y - heights[i]
  r    <- min(bar_r, heights[i] / 2)
  bars <- c(
    bars,
    sprintf(
      paste0(
        '    <path d="M%.2f %.2f V%.2f Q%.2f %.2f %.2f %.2f ',
        'H%.2f Q%.2f %.2f %.2f %.2f V%.2f Z" fill="%s"/>'
      ),
      x, base_y, y + r,
      x, y, x + r, y,
      x + bar_w - r,
      x + bar_w, y, x + bar_w, y + r,
      base_y, fill
    )
  )
}

# ---- assemble ---------------------------------------------------------------

svg <- c(
  sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%d" viewBox="0 0 %.2f %d">',
    W, H, W, H
  ),
  '  <defs>',
  sprintf('    <clipPath id="hex"><polygon points="%s"/></clipPath>', hex_pts),
  '  </defs>',
  sprintf('  <polygon points="%s" fill="%s"/>', hex_pts, col_bg),
  '  <g clip-path="url(#hex)">',
  bars,
  sprintf(
    '    <line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="2.5" stroke-linecap="round" opacity="0.45"/>',
    x0 - 16, base_y + 1.25, x0 + span + 16, base_y + 1.25, col_ink
  ),
  sprintf(
    '    <text x="%.2f" y="424" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="64" font-weight="600" letter-spacing="0.5" fill="%s">skellambrms</text>',
    W / 2, col_ink
  ),
  '  </g>',
  sprintf(
    '  <polygon points="%s" fill="none" stroke="%s" stroke-width="20" clip-path="url(#hex)"/>',
    hex_pts, col_ink
  ),
  '</svg>'
)

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
writeLines(svg, "man/figures/logo.svg")
rsvg::rsvg_png("man/figures/logo.svg", "man/figures/logo.png", width = 520, height = 600)
cat("written\n")
