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

Enzyme.prototype.size=function()
	{
	return this.length();
	}

Enzyme.prototype.at=function(index)
	{
	return this.site.charAt(index);
	}

Enzyme.prototype.toString=function()
	{
	return this.getName();
	}


Enzyme.prototype.getPos5_3=function()
	{
	return this._pos5_3;
	}

Enzyme.prototype.getWeight=function()
	{
	return this._weigth;
	}
	
Enzyme.prototype.isSmall=function()
	{
	return this.getWeight()<5.0;
	}


Enzyme.prototype.getName=function()
	{
	return this.name;
	}
/** the site as declared in rebase, with '^' */
Enzyme.prototype.getDeclaredSite=function()
	{
	return this.raw;
	}

Enzyme.prototype.cut=function(index,plasmid)
	{
	var b=this.at(index);
	switch(plasmid)
		{
		case 'A':return "ANDHVMWR".indexOf(b)!=-1;
                case 'T':return "TNYBDHKW".indexOf(b)!=-1;
                case 'G':return "GNRBDVKS".indexOf(b)!=-1;
                case 'C':return "CNYBHVMS".indexOf(b)!=-1;
                case 'N':return false;
		default: return false;
		}
	}


Enzyme.prototype.getDomId=function()
	{
	return "ENZ"+this.getName();
	}

Enzyme.prototype.getURL=function()
	{
	return "http://rebase.neb.com/rebase/enz/"+this.getName()+".html";
	}

Enzyme.prototype.form=function(div)
	{
	var doc=div.ownerDocument;
	var d1=doc.createElement("div");
	div.appendChild(d1);
	d1.className="enzyme1";
	
	var cbox=doc.createElement("input");
	d1.appendChild(cbox);
	cbox.setAttribute("id",this.getDomId());
	cbox.setAttribute("type","checkbox");
	cbox.setAttribute("checked","true");
	var label=doc.createElement("label");
	label.setAttribute("for",this.getDomId());
	label.appendChild(doc.createTextNode(this.getName()+" ["+this.getDeclaredSite()+"]"));
	d1.appendChild(label);
	
	var a=doc.createElement("a");
	a.setAttribute("title","Link to Rebase");
	a.setAttribute("href",this.getURL());
	a.setAttribute("anchor","_blank");
	a.appendChild(doc.createTextNode("[*]"));
	label.appendChild(a);
	
	var ret={enzyme:null,div:d1,checkbox:cbox};
	ret.enzyme=this;
	return ret;
	}

