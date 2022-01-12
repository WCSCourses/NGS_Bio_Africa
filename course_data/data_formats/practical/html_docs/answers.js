

jQuery.fn.stage_answers = function(options) {

    var hide_answers = 1;
    if ( !hide_answers ) return;

    return this.each(function() {      
        $('li.A').each(function(i, obj) {
            var html = $(this).clone();
            jQuery.data(obj,'ori',html);
            $(this).html('');
            $(this).toggleClass('W');
        });
        $('input').each(function(i,obj) {
            $(this).value = "";
        });
    });
};

function show_answer(event) {
    var obj = event.data.obj;
    $(this).off('click')
    $(obj).empty().append(event.data.ori.contents());
    $(obj).toggleClass('W',false);
    $(obj).toggleClass('A',true);
}

function unlock_answers(obj,section,exp) {
    if ( obj.value!=exp ) 
    {
        alert("Incorrect key");
        return;
    }
    $('li.'+section).each(function(i, obj) {
        var data = jQuery.data(obj,'ori');
        $(this).empty().append(data.contents());
        $(this).toggleClass('W',false);
        $(this).toggleClass('A',true);
    });
}


