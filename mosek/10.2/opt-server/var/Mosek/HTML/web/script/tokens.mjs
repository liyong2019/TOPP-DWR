//-*-javascript-*-
import * as optapi from './optapi.mjs';
import * as utils from './utils.mjs';
import * as gui from './gui.mjs';


function fmtTimeCasual(Tms) {
    var secs_year = 60*60*24*365;
    var secs_day  = 60*60*24;
    var secs_hr   = 60*60;
    var secs_min  = 60;

    var T = Tms < 0 ? -Tms/1000 : Tms/1000;
    var res = [];

    if (T < secs_min)
        return "now";

    if (T > secs_year) {
        var r = Math.floor(T/secs_year);
        res[res.length] = r+(r > 1 ? " years" : " year")
        T = T % secs_year;
    }
    if (T > secs_day) {
        var r = Math.floor(T/secs_day);
        res[res.length] = r+(r > 1 ? " days" : " day")
        T = T % secs_day;
    }
    if (T > secs_hr) {
        var r = Math.floor(T/secs_hr);
        res[res.length] = r+(r > 1 ? " hours" : " hour")
        T = T % secs_hr;
    }
    if (T > secs_min) {
        var r = Math.floor(T/secs_min);
        res[res.length] = r+(r > 1 ? " mins" : " min")
        T = T % secs_min;
    }

    if (res.length == 1) return res[0]+(Tms < 0 ? " ago" : "");
    else return res[0] + ", " + res[1] +(Tms < 0 ? " ago" : "");
}

async function deleteByRow(tr) {
    if (tr.nodeName == "TR") {
        try {
            await optapi.deleteToken(tr.getAttribute("data-serial"));
            tr.remove();
        }
        catch (err) {
        }
    }
}

// async function token_requestCreate(expiry,perms,userid,desc,addrow_func) {
//     if (!isNaN(expiry) && userid.length != 0) {
//         var data = { "Owner"       : userid,
//                      "Expires"     : Math.round(expiry),
//                      "Permissions" : perms,
//                      "Description" : desc };
//         var r = await api_async_createToken(data);
//         var resfield = document.getElementById("output-token");
//         resfield.value = r.Name;
//         addrow_func(r);
//     }
// }


function addRow(tbody,cells,r) {
    var tr = document.createElement("tr");
    tr.setAttribute("data-serial",r.Serial);
    if (r.Expired)
        tr.setAttribute('class','expired');

    for (var j = 0; j < cells.length; ++j) {
        var td  = document.createElement("td");
        switch (cells[j]) {
        case "desc":
	    td.appendChild(document.createTextNode(r.Description));
            break;
        case "owner":
	    td.appendChild(document.createTextNode(r.Owner));
            break;
	case "expiry":
	    td.appendChild(document.createTextNode(fmtTimeCasual(r.Expires*1000-Date.now()), r.Expires));
            td.setAttribute('data-expiry-sec',''+r.Expires);
            break;
        case "perm-submit":
	    td.appendChild(document.createTextNode(r.Permissions['submit'] ? "yes" : "no"));
            break;
        case "perm-admin":
	    td.appendChild(document.createTextNode(r.Permissions['admin'] ? "yes" : "no"));
            break;
	case "ops":
            td.appendChild(
                gui.iconbar_make([
                    gui.icon_make("@remove","Delete Token",function(ev) { deleteByRow(tr); } )
                ]));
            break;
        default:
            break;
        }
	tr.appendChild(td);
    }

    tbody.appendChild(tr);
}

function fmtPermissions(v) {
    let perms = new Array();
    if ((v & 0x0004) != 0) perms.push("admin");
    if ((v & 0x0008) != 0) perms.push("submit");
    if ((v & 0x0010) != 0) perms.push("API");
    if ((v & 0x0001) == 0)
        return "INVALID ("+perms.join(",")+")"
    else
        return perms.join(",");
}
function fmtDateFromUNIX(v) {
    return `${new Date(v).toISOString()}`;
}
export async function initialize(user) {
    let table = await $('#tokens-table').DataTable({
        //dom : 'lB<t>ip',
        // buttons : [
        //     { text : "Delete selected",
        //       action : () => {
        //           let rows = table.rows({ selected : true });
        //           rows.deselect();
        //           Promise.all(rows.data().map(rowdata => optapi.deleteToken(rowdata[0]))).then(() => table.ajax.reload());
        //       }
        //     },
        // ],
        serverSide : true,
        searching  : false,
        paging     : true,
        pageLength : 25,
        columnDefs : [ { targets : 'data-permissions', "render" : (data,type,row,meta) => { return fmtPermissions(data); } },
                       { targets : 'data-expires',     "render" : (data,type,row,meta) => { return fmtDateFromUNIX(data*1000); } },
                       { targets : 'data-ops',         "defaultContent" : "<button>Invalidate</button>" },
                     ],
        select  : false,
        ajax : {
            url         : "/users/api/data",
            type        : "POST",
            contentType : "application/json",
            dataType    : "json",
            data        : (d) => { d['request'] = 'access-tokens';
                                   if (user) d['userid'] = user['Userid'];
                                   return JSON.stringify(d);
                                 }
        }
    });
    $('#tokens-table tbody').on('click','button', async (ev) => {
        console.log($(ev.target).parent().parent());
        let serial = table.row($(ev.target).parent().parent()).data()[0];
        // console.log(table,"---",$(this),$(this).parents('tr'));
        // let serial = table.row( $(this).parents('tr') ).data()[0];
        await optapi.invalidateToken(serial);
        await table.ajax.reload();
    });


    let btn_create = document.getElementById("btn-create-token");
    btn_create.onclick = async function() {
        var expirydate = document.getElementById("input-expiry");
        var expirydays = document.getElementById("input-expiry-days");

        var expirysecs = undefined;
        if (expirydays) {
            var expiry = parseInt(expirydays.value);
            if (isNaN(expiry))
                expirysecs = undefined;
            else
                expirysecs = expiry * 60*60*24;
        } else if (expirydate) {
            expirysecs = Date.parse(expirydate.value);
            if (isNaN(expirysecs))
                expirysecs = undefined;
            else
                expirysecs = Math.round((expirysecs - Date.now()) / 1000);
        }
        let useridnode = document.getElementById("inp-userid");
        let userid = useridnode ? useridnode.value : undefined;


        if (expirysecs !== undefined) {
            var desc   = document.getElementById("inp-desc").value;
            var perms = []
            if (document.getElementById("inp-perm-admin").checked) perms[perms.length]  = "admin";
            if (document.getElementById("inp-perm-submit").checked) perms[perms.length] = "submit";
            if (document.getElementById("inp-perm-api").checked) perms[perms.length]    = "use-api";

            let data = { expires : Math.round(expirysecs),
                         perms   : perms,
                         desc    : desc };
            if (userid !== undefined) data.owner = userid;
            let r = await optapi.createToken(data);
            let resfield = document.getElementById("output-token");
            resfield.value = r.Name;

            table.ajax.reload();
        }
    };
    // var table = document.getElementById("tokens-table");
    // if (table) {
    //     var global = utils.hasClass(table,"global");
    //     let {tbody,thead} = utils.tableElements(table);
    //     let theadrow = thead.firstElementChild

    //     let cells = [];
    //     let statuscell = undefined;
    //     if (theadrow)
    //         for (var n = theadrow.firstElementChild; n ; n = n.nextElementSibling)
    //     	if (n.nodeName == 'TH') {
    //                 var item = n.getAttribute("data-item")
    //                 if (item == 'status')
    //     		statuscell = cells.length;
    //                 cells.push(item);
    //     	}


    //     let data = await optapi.listTokens({user:(global ? undefined : user.Userid)});
    //     while (tbody.firstChild)
    //         tbody.firstChild.remove();
    //     for (var i = 0; i < data.length; ++i)
    //         addRow(tbody,cells,data[i]);

    //     var l = document.getElementsByClassName("active-rows");
    //     for (var i = 0; i < l.length; ++i) {
    //         if (l[i].nodeName == "TABLE")
    //     	gui.table_makeActive(l[i]);
    //     }

    //     document.getElementById("btn-delete-selected").onclick = function() {
    //         utils.table_forallSelected(table,deleteByRow);
    //     };
    // }
}
