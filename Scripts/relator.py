def relator(edge):
  a = edge.get_arrow()
  start = a.copy()
  relator = []
  while 1:
    relator.append(a.glued())
    a.next()
    if a == start:
      break
  return relator