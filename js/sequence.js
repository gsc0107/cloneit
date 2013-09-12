function Sequence()
	{
	this.name="";
	this.sequence="";
	}

Sequence.prototype.getName=function()
	{
	return this.name;
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

Sequence.prototype.toString=function(index)
	{
	return this.getName()+" length:"+this.size();
	}

Sequence.create=function(name,sequence)
	{
	var S=new Sequence();
	S.name=name;
	S.sequence=sequence.replace(/\s+/g,'').toUpperCase();
	
	return S;
	}
