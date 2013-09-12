function Rebase()
	{
	this.enzymes=[];
	}

Rebase.prototype.size=function()
	{
	return this.enzymes.length;
	}

Rebase.prototype.toString=function()
	{
	return "Rebase N="+this.size();
	}

Rebase.prototype.get=function(index)
	{
	if(index<0 || index >= this.size()) throw "Rebase: index out of range 0<="+index+"<"+this.size();
	return this.enzymes[index];
	}

Rebase.prototype.add=function(name,size,weight)
	{
	var enz=new Enzyme(name,size,weight);
	this.enzymes.push(enz);
	return enz;
	}

Rebase.uncheckSmalls=function(elementsForms,weight)
	{
	for(var i=0;i< elementsForms.length; ++i)
		{
		var f=elementsForms[i];
		if(f.checkbox.checked && f.enzyme.getWeight()< weight)
			{
			f.checkbox.checked=false;
			}
		}
	}

Rebase.prototype.form=function(div)
	{
	var i,e;
	var elementsForms=[];
	var doc=div.ownerDocument;

	var but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Unckeck All"));
	but.addEventListener("click",function()
		{
		for(var i=0;i< elementsForms.length; ++i)
			{
			elementsForms[i].checkbox.checked=false;
			}
		});
		
	but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Ckeck All"));
	but.addEventListener("click",function()
		{
		for(var i=0;i< elementsForms.length; ++i)
			{
			elementsForms[i].checkbox.checked=true;
			}
		});
	
	/** 4 **/
	but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Discard < 4"));
	but.addEventListener("click",function()
			{
			Rebase.uncheckSmalls(elementsForms,4);
			});
	/** 5 **/
	but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Discard < 5"));
	but.addEventListener("click",function()
			{
			Rebase.uncheckSmalls(elementsForms,5);
			});
	/** 6 **/
	but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Discard < 6"));
	but.addEventListener("click",function()
			{
			Rebase.uncheckSmalls(elementsForms,6);
			});	
	/** 7 **/
	but=doc.createElement("button");
	div.appendChild(but);
	but.appendChild(doc.createTextNode("Discard < 7"));
	but.addEventListener("click",function()
			{
			Rebase.uncheckSmalls(elementsForms,7);
			});
	
	var d1=doc.createElement("div");
	d1.className="rebase1";
	div.appendChild(d1);
	for(i=0;i< this.size();++i)
		{
		var form=this.get(i).form(d1)
		elementsForms.push(form);
		}
	return elementsForms;
	}
