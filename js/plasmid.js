function Plasmid(sequence,is_forward)
	{
	this.sequence=sequence;
	this.is_forward=is_forward;
	this.sites=[];
	}

Plasmid.prototype.length=function()
	{
	return this.sequence.length();
	}

Plasmid.prototype.size=function()
	{
	return this.length();
	}

Plasmid.prototype.at=function(index)
	{
	index=index%this.size();
	if(!this.is_forward)
		{
		return this.sequence.at(index);
		}
	
	switch(this.sequence.at((this.size()-1)-index))
		{
		case "A": return "T";
		case "T": return "A";
		case "G": return "C";
		case "C": return "G";
		default: return "N";
		}
		
	}
	
Plasmid.prototype.digest=function(rebase)
	{
	var i,k,site;
	this.sites=[];
	
	for(k=0;k< rebase.size();++k)
		{
		
		var j;
		var enz=rebase.get(k);

		var enz_len=enz.size();

		
		for(i=0;i< this.size();++i)
			{
			
			for(j=0;j< enz_len ;++j)
				{
				if(!enz.cut(j,this.at(i+j))) break;
				}

			if(j==enz_len)
				{
				
				site=new Site(this,enz,i);
				this.sites.push(site);
				
				}
			
			}
		}
	this.sites.sort(Site.sortOnPos);
	}

Plasmid.prototype._lower_bound=function(beg,end,pos)
	{
	var half,middle,len;
	len = last - first;
	while (len > 0)
		    {
		    half = len / 2;
		    middle = first + half;
		    
		    if ( this.sites[middle].pos < 0  )
		            {
		            first = middle + 1;
		            len -= half + 1;
		            }
		    else
		            {
		            len = half;
		            }
		    }
	return first;
	}

