(window.webpackJsonp=window.webpackJsonp||[]).push([[0],{131:function(e,t,a){"use strict";a.r(t);var n=a(0),r=a.n(n),l=a(41),i=a.n(l),o=(a(54),a(132)),c=a(139),m=a(4),s=a(25),u=a(140),d=Object(s.a)(function(e){var t=e.theme,a=Object(m.a)(e,["theme"]);return r.a.createElement(c.a,Object.assign({display:"flex",padding:16,background:t.palette.blue.base},a),r.a.createElement(c.a,{flex:1,alignItems:"center",display:"flex"},r.a.createElement(u.a,{color:t.colors.border.default},"Capgen")),r.a.createElement(c.a,null))}),p=a(141),f=a(142),E=a(16),b=a(30),g=a.n(b),h=a(44),x=a(3),v=a(43),y=a(5),O=a(7),j=a(6),w=a(8),k=window.libcapgen(),R=null,z=[];k.onRuntimeInitialized=function(){R={optimize:function(){return k.optimize.apply(k,arguments)}},z.map(function(e){return e()})};var B=function(){return new Promise(function(e){null===R?z.push(function(){return e(R)}):e(R)})},F=a(146),S=function(e){function t(){var e,a;Object(y.a)(this,t);for(var n=arguments.length,r=new Array(n),l=0;l<n;l++)r[l]=arguments[l];return(a=Object(O.a)(this,(e=Object(j.a)(t)).call.apply(e,[this].concat(r)))).state={history:[],loading:!1,metrics:{totalReturn:0,risk:0}},a.optimize=Object(v.a)(g.a.mark(function e(){var t,n,r,l=arguments;return g.a.wrap(function(e){for(;;)switch(e.prev=e.next){case 0:return e.next=2,B();case 2:t=e.sent,n=t.optimize.apply(t,l),r={totalReturn:n.totalReturn,risk:n.risk},a.setState(Object(x.a)({},a.state,{history:[{title:"Backtest",timestamp:Object(F.a)(new Date,"yyyy-MM-dd HH:mm:ss"),metrics:Object(x.a)({},r)}].concat(Object(h.a)(a.state.history)),metrics:Object(x.a)({},r)}));case 6:case"end":return e.stop()}},e)})),a}return Object(w.a)(t,e),t}(E.a),P=Object(s.a)(function(e){var t=e.theme,a=Object(m.a)(e,["theme"]);return r.a.createElement(c.a,Object.assign({height:"100%",overflowY:"auto",width:"200px",background:"tint1",borderRight:"1px solid ".concat(t.colors.border.default)},a),r.a.createElement(E.c,{to:[S]},function(e){return e.state.history.map(function(e,a){return r.a.createElement(c.a,{key:a,padding:16,borderBottom:"1px solid ".concat(t.colors.border.default)},r.a.createElement(u.a,{size:200},e.title),r.a.createElement(u.a,{size:100},e.timestamp),r.a.createElement(p.a,{marginBottom:0},r.a.createElement(X,{icon:"symbol-triangle-up"},e.metrics.totalReturn.toFixed(2)),r.a.createElement(X,{icon:"symbol-triangle-down"},e.metrics.risk.toFixed(2))))})}))}),X=function(e){var t=e.children,a=Object(m.a)(e,["children"]);return r.a.createElement(f.a,Object.assign({size:100,display:"inline-block",marginTop:0,marginBottom:0,marginRight:32,paddingLeft:0,color:"muted",fontSize:11},a),t)},W=a(145),C=a(24),D=a(147),I=a(143),Y=a(144),M=function(e){var t=e.children,a=Object(m.a)(e,["children"]);return r.a.createElement(c.a,Object.assign({marginTop:16},a),t)},A=function(e){var t=e.children,a=Object(m.a)(e,["children"]);return r.a.createElement(I.a,a,t)},T=function(e){var t=e.children,a=Object(n.useState)(0),l=Object(C.a)(a,2),i=l[0],o=l[1];return r.a.createElement(r.a.Fragment,null,r.a.createElement(Y.a,{marginX:-4},t.map(function(e,t){var a=e.props,n=a.title,l=(a.children,Object(m.a)(a,["title","children"]));return r.a.createElement(A,Object.assign({key:t,is:"a",isSelected:t===i,onSelect:function(){return o(t)}},l),n)})),t.map(function(e,t){return i===t&&r.a.createElement(M,{key:t},e.props.children)}))},H={populationSize:100,elitismCopies:2,generations:100,mutationRate:.02,crossoverRate:.02,steps:4,riskAversion:.5,initialFundingRatio:1.3,targetFundingRatio:1.3},J=Object(s.a)(function(e){var t=e.theme,a=Object(n.useState)(Object(x.a)({},H)),l=Object(C.a)(a,2);l[0],l[1];return r.a.createElement(r.a.Fragment,null,r.a.createElement(c.a,{background:"tint1",borderBottom:"1px solid ".concat(t.colors.border.default),display:"flex",paddingX:16,paddingY:8,flexDirection:"column"},r.a.createElement(c.a,{flex:1},r.a.createElement(u.a,{size:200},"Parameters"))),r.a.createElement(c.a,{borderBottom:"1px solid ".concat(t.colors.border.default),display:"flex",padding:16,flexDirection:"column"},r.a.createElement(T,null,r.a.createElement(A,{title:"Portfolio"},r.a.createElement(G,null)),r.a.createElement(A,{title:"Simulation"},r.a.createElement(K,null)),r.a.createElement(A,{title:"Optimizer"},r.a.createElement(L,null)),r.a.createElement(A,{title:"Correlations",disabled:!0},"Correlations"))))}),N=function(e){return r.a.createElement(D.a,Object.assign({flexBasis:"28%",marginX:8},e))},G=function(){return r.a.createElement(c.a,{display:"flex",flexWrap:"wrap",marginX:-8},r.a.createElement(N,{label:"Risk aversion",hint:"Set to 0.0 to optimize without considering risk.",placeholder:"0.5",value:.5}),r.a.createElement(N,{label:"Initial funding ratio",hint:"Initial funding ratio of the portfolio.",placeholder:"1.3",value:1.3}),r.a.createElement(N,{label:"Target funding ratio",hint:"Penalize portfolios with funding ratio below this value.",placeholder:"1.3",value:1.3}))},K=function(){return r.a.createElement(c.a,{display:"flex",flexWrap:"wrap",marginX:-8},r.a.createElement(N,{label:"Steps",hint:"Number of steps to simulate.",placeholder:"4",value:4}))},L=function(){return r.a.createElement(c.a,{display:"flex",flexWrap:"wrap",marginX:-8},r.a.createElement(N,{label:"Population size",hint:"Number of individuals in the population.",placeholder:"100",value:100}),r.a.createElement(N,{label:"Elitism copies",hint:"Keep a number of clones of the best individuals at iteration.",placeholder:"2",value:2}),r.a.createElement(N,{label:"Generations",hint:"Iterate through this many generations.",placeholder:"100",value:100}),r.a.createElement(N,{label:"Mutation rate",hint:"Probability of mutating a gene.",placeholder:"0.02",value:.02}),r.a.createElement(N,{label:"Crossover rate",hint:"Probability of crossing two individuals.",placeholder:"0.02",value:.02}))},$=Object(s.a)(function(e){var t=e.theme;return r.a.createElement(r.a.Fragment,null,r.a.createElement(c.a,{background:"tint1",borderBottom:"1px solid ".concat(t.colors.border.default),display:"flex",paddingX:16,paddingY:8,flexDirection:"column"},r.a.createElement(c.a,{flex:1},r.a.createElement(u.a,{size:200},"Results"))),r.a.createElement(c.a,{display:"flex",padding:16,flexDirection:"column"},r.a.createElement(T,null,r.a.createElement(A,{title:"Metrics"},r.a.createElement(c.a,{display:"flex",flexWrap:"wrap",marginX:-8,marginY:-8},r.a.createElement(E.c,{to:[S]},function(e){return r.a.createElement(r.a.Fragment,null,r.a.createElement(q,{title:"Expected Return",value:e.state.metrics.totalReturn.toFixed(2)}),r.a.createElement(q,{title:"Expected Risk",value:e.state.metrics.risk.toFixed(2)}))}))),r.a.createElement(A,{title:"Performance"},"Performance"),r.a.createElement(A,{title:"Weights"},"Weights"))))}),q=function(e){var t=e.title,a=e.value,n=e.children,l=Object(m.a)(e,["title","value","children"]);return r.a.createElement(c.a,Object.assign({background:"tint1",padding:16,flexBasis:"28%",marginX:8,marginY:8,border:"default"},l),r.a.createElement(u.a,{size:200,marginBottom:8},t),"undefined"===typeof a?n:r.a.createElement(u.a,{size:800},a))},Q=Object(s.a)(function(e){var t=e.theme,a=Object(m.a)(e,["theme"]);return r.a.createElement(c.a,Object.assign({display:"flex",flexDirection:"column"},a),r.a.createElement(c.a,{borderBottom:"1px solid ".concat(t.colors.border.default),display:"flex",padding:16,alignItems:"center"},r.a.createElement(c.a,{flex:1},r.a.createElement(u.a,null,"Backtest")),r.a.createElement(c.a,null,r.a.createElement(E.c,{to:[S]},function(e){return r.a.createElement(W.a,{appearance:"primary",iconBefore:"play",onClick:function(){return e.optimize({populationSize:10,elitismCopies:2,generations:10,mutationRate:.02,crossoverRate:.02,steps:4,riskAversion:.5,initialFundingRatio:1.3,targetFundingRatio:1.3})}},"Run")}))),r.a.createElement(J,null),r.a.createElement($,null))}),U=function(){return r.a.createElement(E.c,{to:[S]},function(e){return e.loading?r.a.createElement(o.a,{transform:"translate(-50%, -50%)",position:"absolute",left:"50%",top:"50%"}):r.a.createElement(c.a,{height:"100vh",display:"flex",flexDirection:"column"},r.a.createElement(d,null),r.a.createElement(c.a,{display:"flex",height:"calc(100% - ".concat("52px",")")},r.a.createElement(c.a,{flex:0,height:"100%"},r.a.createElement(P,null)),r.a.createElement(c.a,{flex:1},r.a.createElement(Q,null))))})};Boolean("localhost"===window.location.hostname||"[::1]"===window.location.hostname||window.location.hostname.match(/^127(?:\.(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)){3}$/));i.a.render(r.a.createElement(E.b,null,r.a.createElement(U,null)),document.getElementById("root")),"serviceWorker"in navigator&&navigator.serviceWorker.ready.then(function(e){e.unregister()})},48:function(e,t,a){e.exports=a(131)},54:function(e,t,a){}},[[48,1,2]]]);
//# sourceMappingURL=main.8f0077c4.chunk.js.map