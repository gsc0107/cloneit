function Site(sequence,enzyme,pos)
	{
	this.sequence=sequence;
	this.enzyme=enzyme;
	this.pos=pos;
	}


Site.prototype.getSequence=function()
	{
	return this.sequence;
	}

Site.prototype.getEnzyme=function()
	{
	return this.enzyme;
	}

Site.prototype.getPosition=function()
	{
	return this.pos;
	}

Site.prototype.length=function()
	{
	return this.getEnzyme().length();
	}


Site.prototype.at=function(index)
	{
	return this.getSequence().charAt(this.getPosition()+index);
	}

Site.compatible=function(L,R,config)
	{
	return true;
	}


Site.sortOnPos=function(a,b)
	{
	return a.getPosition() - b.getPosition();
	}


