function grab(id)
{
  return document.getElementById(id);
}



function set_value(id, val)
{
  grab(id).value = val;
}



function get_float(id, def)
{
  var x = parseFloat( grab(id).value );
  return !isNaN(x) && isFinite(x) ? x : def;
}



function get_int(id, def)
{
  var x = parseInt( grab(id).value );
  return !isNaN(x) && isFinite(x) ? x : def;
}



/* check if `n' is a valid integer */
function is_int(n)
{
  return !isNaN(parseInt(n)) && isFinite(n);
}



/* check if `n' is a valid floating number */
function is_float(n)
{
  return !isNaN(parseFloat(n)) && isFinite(n);
}



