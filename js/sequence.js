function Sequence()
	{
	this.name="";
	this.sequence="";
	}

Sequence.prototype.length=function()
	{
	return this.sequence.length;
	}

Sequence.prototype.size=function()
	{
	return this.length();
	}

Sequence.prototype.at=function(index)
	{
	return this.sequence.charAt(index);
	}
