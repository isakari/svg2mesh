{
  # gray(整数) → gray(小数)
  if ($0 ~ /gray\([0-9]+\)/ && $0 !~ /\./) {
    g = $0
    sub(/^.*gray\(/, "", g)
    sub(/\).*$/, "", g)
    val = g + 0
    perc = sprintf("%.4f", val / 255 * 100)
    sub(/gray\([0-9]+\)/, "gray(" perc ")")
  }

  # gray(小数%) → gray(小数)
  if ($0 ~ /gray\([0-9]+\.[0-9]+%\)/) {
    g = $0
    sub(/^.*gray\(/, "", g)
    sub(/%\).*$/, "", g)
    sub(/gray\([0-9]+\.[0-9]+%\)/, "gray(" g ")")
  }

  print
}
