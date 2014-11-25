
/* print array */
function parr(a, name)
{
  var n = a.length;
  var b = new Array(n); // copy the array to avoid the delay in Firefox/Chrome
  cparr(b, a, n);
  console.log(name, b);
}



/* print vector */
function pvec(a, name)
{
  parr(a, name);
}



/* print matrix */
function pmat(a, name)
{
  for ( var j = 0; j < a.length; j++ )
    parr(a[j], name + "[" + j + "]");
}
