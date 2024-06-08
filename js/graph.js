// https://gist.github.com/heybignick/3faf257bbbbc7743bb72310d03b86ee8

var color = d3.scaleOrdinal(d3.schemeCategory20);
var radius = 5;

var simulation = d3.forceSimulation()
    .force("link", d3.forceLink()
      .id(function(d) { return d.id; })
      .distance(function(d) { return d.distance; })
      )
    .force("charge", d3.forceManyBody()
      .strength(-30)
      .distanceMax(100))
    .force("center", d3.forceCenter(width / 2, height / 2));

r2d3.onRender(function(graph, svg, width, height, options) {
  var link = svg.append("g")
    .attr("class", "links")
    .selectAll("line")
    .data(graph.edges)
    .enter().append("line")
      .attr("stroke-width", 1)
      .attr("stroke", "#999")
      .attr("stroke-opacity", 0.6);

  var node = svg.append("g")
      .attr("class", "nodes")
      .selectAll("g")
      .data(graph.nodes)
      .enter().append("g");
    
  var circles = node.append("circle")
    .attr("r", radius)
    .attr("fill", function(d) { return color(d.group); });
  
  var drag_handler = d3.drag()
    .on("start", dragstarted)
    .on("drag", dragged)
    .on("end", dragended);
    
  drag_handler(node);
  
  var labels = node.append("text")
      .text(function(d) { return d.label; })
      .attr('x', 6)
      .attr('y', 3)
      .attr("font-family", "sans-serif")
      .attr("font-size", "10px");
      
  node.append("title")
      .text(function(d) { return d.id; });

  simulation
      .nodes(graph.nodes)
      .on("tick", ticked);

  simulation.force("link")
      .links(graph.edges);

  function ticked() {
    //constrains the nodes to be within a box
    /*
    node.attr("cx", function(d) { 
      return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
        .attr("cy", function(d) { 
          return d.y = Math.max(radius, Math.min(height - radius, d.y)); });
    */
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    //node.attr("cx", function(d) { return d.x; })
    //    .attr("cy", function(d) { return d.y; });
    node.attr("transform", function(d) {
      return "translate(" + d.x + "," + d.y + ")";
    })
  }
});

function dragstarted(d) {
  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
  d.fx = d.x;
  d.fy = d.y;
}

function dragged(d) {
  d.fx = d3.event.x;
  d.fy = d3.event.y;
}

function dragended(d) {
  if (!d3.event.active) simulation.alphaTarget(0);
  d.fx = null;
  d.fy = null;
}

