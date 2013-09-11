function Enzyme(name,site,commercial)
	{
	var i=0;
	this.name=name;
	this.raw=site;
	this.commercial=commercial;
	this.site="";
	this._pos5_3=-1;
	this._weigth=0.0;
	for(i=0;i< site.length;++i)
		{
		var c=site.charAt(i);
		if(c=='^')
			{
			this._pos5_3=this.site.length;
			}
		else
			{
			switch(c)
				{
				case "A": case "T": case "G": case "C": this._weigth+=1.0; break;
				case "S": case "W": case "Y": case "R":  case "M": case "K":  this._weigth+=0.5; break;
				case "B": case "D": case "H": case "V": this._weigth+=0.75; break;
				case "N": break;
				default:break;
				}
			
			this.site+=c;
			}
		}
	}

Enzyme.prototype.length=function()
	{
	return this.site.length;
	}

Enzyme.prototype.at=function(index)
	{
	return this.site.charAt(index);
	}

Enzyme.prototype.getPos5_3=function()
	{
	return this._pos5_3;
	}

Enzyme.prototype.getWeight=function()
	{
	return this._weigth;
	}

Enzyme.prototype.getName=function()
	{
	return this.name;
	}


Enzyme.prototype.form=function(doc)
	{
	var df=doc.createDocumentFragment();
	
	var cbox=doc.createElement("input");
	cbox.setAttribute("type","checkbox");
	var label=doc.createElement("label");
	label.appendChild(doc.createTextNode(this.getName()));
	
	return df;
	}

