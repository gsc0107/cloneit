var App={
	rebase:null,
	enzymeCheckboxes:null,
	vector:null,
	plasmid:null
	};

this.addEventListener("load",function(evt){
 	var r=document.getElementById("rebase");
 	App.rebase = makeRebase();
 	App.enzymeCheckboxes = App.rebase.form(r);
 	
 	r=document.getElementById("x1");
 	r.value=pgbt9.sequence.replace(/(.{50})/g, "$1\n");
 	
 	App.vector=pgbt9;
 	App.plasmid=new Plasmid(App.vector,true);
 	
 	r=document.getElementById("run");
 	
 	r.addEventListener("click",function()
		{
		
		App.plasmid.digest(App.rebase);
		r=document.getElementById("result");
		if(App.plasmid.sites.length==0) return;
		var i=0;
		for(i=0;i< App.plasmid.sites.length;i++)
			{
			var enz1=App.plasmid.sites[i];
			var enz2=App.plasmid.sites[i+1==App.plasmid.sites.length?0:i+1];
			
			var tr=document.createElement("tr");
			r.appendChild(tr);
			
			var td=document.createElement("td");tr.appendChild(td);
			td.appendChild(document.createTextNode(enz1.getEnzyme().getName()));
			
			td=document.createElement("td");tr.appendChild(td);
			td.appendChild(document.createTextNode(enz1.getPosition()));
			
			td=document.createElement("td");tr.appendChild(td);
			td.appendChild(document.createTextNode(enz2.getPosition()));
			td=document.createElement("td");tr.appendChild(td);
			td.appendChild(document.createTextNode(enz2.getEnzyme().getName()));
			}
		alert("ok");
		});
	
 	},false);	
