function check_change() {
  
  /* start from root */
  var queue = [TREE['root']];
  while (queue) {
    var current = queue.pop();
    
    // get the children
    var children = TREE[current]['children'];

    // check for current char
    var this_node = document.getElementById(TREE[current]['label']+':label');
    var this_content = this_node.dataset.value;

    // iterate over children
    for (var i=0,child; child=children[i]; i++) {
      // get child value
      var child_node = document.getElementById(TREE[child]['label']+':label');
      var child_content = child_node.dataset.value;
      
      // change border if child not equals mother
      // change border if child not equals mother
      if (child_content != this_content) {
	if (TREE[child]['leave']) {
	  child_node.childNodes[0].style.border = '10px double black';
	}
	else {
	  child_node.style.border = '10px double black';
	}
      }

      queue.push(child);
    }
  }
}



