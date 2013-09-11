function Plasmid(sequence,is_forward)
	{
	this.sequence=sequence;
	this.is_forward=is_forward;
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
	var i,j,k,enz,site;
	this.sites=[];
	for(k=0;k< rebase.size();++k)
		{
		enz=rebase.get(k);
		for(i=0;i< this.size();++i)
			{
			for(j=0;j<enz.size();++j)
				{
				
				}
			if(j==enz.size())
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

