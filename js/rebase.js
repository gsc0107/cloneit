function Rebase()
	{
	this.enzymes=[];
	}

Rebase.prototype.size=function()
	{
	return this.enzymes.length;
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

Rebase.prototype.form=function(doc)
	{
	var i,e;
	var d1=document.createElement("div");
	for(i=0;i< this.size();++i)
		{
		var e = document.createElement("");
		}
	}
