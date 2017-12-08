/**
    @Author Yinchen Xu yxu72@illinois.edu
    Reference: UIUC CS 519 Fall 2017 source codes, by Prof. Eric Shaffer
*/

//Math functions

//Simple gaussian function, return a value based on the canvas coordinate position.
function gaussian(pt){
	return Math.exp(-(pt[0] * pt[0] + pt[1] * pt[1]));
}

//Gaussian function's gradient, return a vector based on the canvas coordinate position
function gaussian_gradient(pt){
  var dx = -2*pt[0]*gaussian(pt);
  var dy = -2*pt[1]*gaussian(pt);
	return [dx,dy];
}

//Helper function to normalize a 2D vector.
function normalize2D(v){
  var len = Math.sqrt(v[0] * v[0] + v[1] * v[1]);
  if (len == 0.0)
    {
       console.log("Zero length gradient");
       return ([0.0,0.0]);
    }
  return [v[0]/len,v[1]/len];
}

/**
    Helper function to do euler integration
    @param pt - canvas coordinate position of the input
    @param h - step size
    @param steps - total steps of this euler integration
    @param get_vector - function to get the corresponding vector in the vector field
*/
function euler_integration(pt, h, radius, get_vector)
{
    var result_numerator = 0.0;
    var result_denominator = 0.0;
    var ln = [[pt[0], pt[1]]];
    var last = [pt[0], pt[1]];

    for(i=0;radius > 0 && inBound(last) && i < 450;i++)
    {
        v = get_vector(ln[i]);
        ln.push( [ln[i][0]+h*v[0],ln[i][1]+h*v[1]]);

        var pixel_pt = pt2pixel(400, 400, [-1., 1.], [-1., 1.], ln[i][0], ln[i][1]);

        result_denominator += Math.exp(-(blur_radius - radius)*(blur_radius - radius));
        result_numerator += Math.exp(-(blur_radius - radius)*(blur_radius - radius)) * noise_greyscales[pixel_pt[0]][pixel_pt[1]];

        radius -= step_len * Math.sqrt(v[0]*v[0] + v[1] * v[1]);
        last = [ln[i][0]+h*v[0],ln[i][1]+h*v[1]];

    }

    return [result_numerator, result_denominator];
}

function euler_integration_reverse(pt, h, radius, get_vector)
{
    var result_numerator = 0.0;
    var result_denominator = 0.0;
    var ln = [[pt[0], pt[1]]];
    var pixel_pts = [pt2pixel(400, 400, [-1., 1.], [-1., 1.], pt[0], pt[1])];
    var last = [pt[0], pt[1]];
    for(i=0; radius > 0 && inBound(last) && i < 450; i++)
    {
        v = get_vector(ln[i]);
        ln.push( [ln[i][0]-h*v[0],ln[i][1]-h*v[1]] );
        var pixel_pt = pt2pixel(400, 400, [-1., 1.], [-1., 1.], ln[i][0], ln[i][1]);

        if(i!=0 && inBound_pixel(pixel_pt))
        {
            result_denominator += Math.exp(-(blur_radius - radius)*(blur_radius - radius));
            result_numerator += Math.exp(-(blur_radius - radius)*(blur_radius - radius)) * noise_greyscales[pixel_pt[0]][pixel_pt[1]];
        }
        radius -= step_len * Math.sqrt(v[0]*v[0] + v[1] * v[1]);
        last = [ln[i][0]-h*v[0],ln[i][1]-h*v[1]];

    }

    return [result_numerator, result_denominator];
}

//render functions and variables

var x_extent = [-1.0, 1.0];
var y_extent = [-1.0, 1.0];
var red_color = [255.0, 0.0, 0.0, 255.0];
var noise_colors = [];
var vector_func = gaussian_gradient;
var blur_radius = 0.05;
var hedge_res = 10;
var hedge_scalar = 1.0;

var step_lens = [0.0002, 0.0004];//Can be modified to balance performance and accuracy, lens[0] for middle area(small vector values), lens[1] for surrounding area (large vector values);

function main()
{
    render();

}

function render()
{
    console.log(vector_func([0.0, 0.0]));
    console.log(vector_func([0.1, 0.1]));

    generate_noise();
    generate_LIC();
    blur_radius = parseFloat(document.getElementById("blur_radius").value);
    hedge_res = parseFloat(document.getElementById("hedge_res").value);
    hedge_scalar = parseFloat(document.getElementById("hedge_scale").value);
    if (document.getElementById("show_hedge").checked)
        draw_hedge(); 
}

var sample_locations_vec = [];

function inBound(pt)
{
    return pt[0] >= x_extent[0] && pt[0] <= x_extent[1] && pt[1] >= y_extent[0] && pt[1] <= y_extent[1];
}

function inBound_pixel(pixel_pt)
{
    return pixel_pt[0] > 0 && pixel_pt[0] < 400 && pixel_pt[1] > 0 && pixel_pt[1] < 400;
}

function draw_hedge()
{
    if(!document.getElementById("random_hedge").checked)
    {
        sample_locations_vec = [];
        for(var i = 0; i < hedge_res; i++)
        {
            sample_locations_vec.push([]);
            for(var j = 0; j < hedge_res; j++)
            {
                sample_locations_vec[i].push(vector_func([-1.0 + 2 * (i+1) / (hedge_res+1), -1.0 + 2 * (j+1)/(hedge_res + 1)]));
            }
        }

        for(var i = 0; i < hedge_res; i++)
        {
            for(var j = 0; j < hedge_res; j++)
            {
                draw_hedge_line(i,j,sample_locations_vec[i][j]);
            }
        }
    }
    else
    {
        var random_locations = [];
        for(var i = 0; i < hedge_res; i++)
        {
            random_locations.push([]);
            for(var j = 0; j < hedge_res; j++)
            {
                random_locations[i].push([-1. + 2 * Math.random(), -1. + 2 * Math.random()]);
            }
        }
        sample_locations_vec = [];
        for(var i = 0; i < hedge_res; i++)
        {
            sample_locations_vec.push([]);
            for(var j = 0; j < hedge_res; j++)
            {
                sample_locations_vec[i].push(vector_func(random_locations[i][j]));
            }
        }
        for(var i = 0; i < hedge_res; i++)
        {
            for(var j = 0; j < hedge_res; j++)
            {
                draw_hedge_line_random(random_locations[i][j],sample_locations_vec[i][j]);
            }
        }
    }
}

function draw_hedge_line_random(position, vector)
{
    var canvas = document.getElementById('example');
    var ctx = canvas.getContext('2d');
    var pixel_pt = pt2pixel(400., 400., [-1., 1.], [-1., 1.], position[0], position[1]);
    ctx.beginPath();
    ctx.strokeStyle = "blue";
    ctx.moveTo(pixel_pt[0], pixel_pt[1]);
    ctx.lineTo(pixel_pt[0] + 10 * hedge_scalar * vector[0], pixel_pt[1] + 10 * hedge_scalar * vector[1]);
    ctx.stroke();
}

function draw_hedge_line(row, column, vector)
{
    var canvas = document.getElementById('example');
    var ctx = canvas.getContext('2d');
    var pixel_pt = pt2pixel(400., 400., [-1., 1.], [-1., 1.], -1 + 2*(row+1)/(hedge_res+1), -1 + 2*(column+1)/(hedge_res+1));
    ctx.beginPath();
    ctx.strokeStyle = "red";
    ctx.moveTo(pixel_pt[0], pixel_pt[1]);
    ctx.lineTo(pixel_pt[0] + 10 * hedge_scalar * vector[0], pixel_pt[1] + 10 * hedge_scalar * vector[1]);
    ctx.stroke();
}

function pt2pixel(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];

	var xlen = (p_x-x_extent[0])/(x_extent[1]-x_extent[0]);
  var ylen = (p_y-y_extent[0])/(y_extent[1]-y_extent[0]);

	pt[0]=Math.round(xlen*width);
	pt[1]=Math.round(ylen*height);
	return pt;
}
var img_colors = [];
var counter = 0;
var print_counter = true;
var step_len = step_lens[0];
function generate_LIC()
{
    img_colors = [];
    var canvas = document.getElementById("example");
    var ctx = canvas.getContext('2d');
    var imgData = ctx.getImageData(0,0,canvas.width, canvas.height);
    for(var i = 0; i < canvas.width; i++)
    {
        img_colors.push([]);
        for(var j = 0; j < canvas.height; j++)
        {
            var greyscale_denominator = 0.1;
            var greyscale_numerator = 0.0;
            var canvas_pt = pixel2pt(400., 400., [-1., 1.], [-1., 1.], i, j);
            if((i < 180 || i > 220) && (j < 180 || j > 220))
            {
                step_len = step_lens[1];    
            }
            else
                step_len = step_lens[0];
            var pair = euler_integration(canvas_pt,step_len,blur_radius,gaussian_gradient);
            greyscale_numerator += pair[0];
            greyscale_denominator += pair[1];
            pair = euler_integration_reverse(canvas_pt,step_len,blur_radius,gaussian_gradient);
            greyscale_numerator += pair[0];
            greyscale_denominator += pair[1];
            var greyscale = 0.0;

            if(!document.getElementById("show_color").checked)
            {
                if(greyscale_denominator != 0)
                {
                    greyscale = greyscale_numerator / greyscale_denominator
                    img_colors[i].push(greyscale);
                }
                else
                {
                    img_colors[i].push(0.0);
                }

                var x = (j*canvas.width + i)*4;
                imgData.data[x]=greyscale;
                imgData.data[x+1]= greyscale;
                imgData.data[x+2]= greyscale;
                imgData.data[x+3]= 255.0;
            }
            else
            {
                var curr_vec = gaussian_gradient(canvas_pt);
                var magnitude = Math.sqrt(curr_vec[0] * curr_vec[0] + curr_vec[1] * curr_vec[1]);
                var curr_color = rainbow_colormap(magnitude, 0.0, 0.8578);//max abs(grad(e^(-r^2))) = sqrt(2/e) = 0.8578
                curr_color[0] = curr_color[0] / 255.;
                curr_color[1] = curr_color[1] / 255.;
                curr_color[2] = curr_color[2] / 255.;//The ratio of all color components * greyscale = final color;
                
                if(greyscale_denominator != 0)
                {
                    greyscale = greyscale_numerator / greyscale_denominator
                    img_colors[i].push(greyscale);
                }
                else
                {
                    img_colors[i].push(0.0);
                }

                var x = (j*canvas.width + i)*4;
                imgData.data[x]=greyscale * curr_color[0];
                imgData.data[x+1]= greyscale * curr_color[1];
                imgData.data[x+2]= greyscale * curr_color[2];
                imgData.data[x+3]= 255.0;
                
            }

        }
    }
    ctx.putImageData(imgData, 0, 0);
}

function rainbow_colormap(fval,fmin,fmax){
	var dx=0.8;
	var fval_nrm = (fval-fmin)/(fmax-fmin);
	var g = (6.0-2.0*dx)*fval_nrm +dx;
	var R = Math.max(0.0,(3.0-Math.abs(g-4.0)-Math.abs(g-5.0))/2.0 )*255;
	var G = Math.max(0.0,(4.0-Math.abs(g-2.0)-Math.abs(g-4.0))/2.0 )*255;
	var B = Math.max(0.0,(3.0-Math.abs(g-1.0)-Math.abs(g-2.0))/2.0 )*255;
	color = [Math.round(R),Math.round(G),Math.round(B),255];
	return color;
}


function pixel2pt(width,height,x_extent,y_extent, p_x,p_y){
	var pt = [0,0];
	xlen=x_extent[1]-x_extent[0]
	ylen=y_extent[1]-y_extent[0]
	pt[0]=(p_x/width)*xlen + x_extent[0];
	pt[1]=(p_y/height)*ylen + y_extent[0];
	return pt;
}


//Helper function to generate a noise graph.
var noise_greyscales = []
function generate_noise()
{
    noise_greyscales = [];
    noise_colors = [];
    var canvas_noise = document.getElementById("noise");
    var ctx_noise = canvas_noise.getContext('2d');
    var noiseData = ctx_noise.getImageData(0,0,canvas_noise.width, canvas_noise.height);
    for(var i = 0; i < canvas_noise.width; i++)
    {
        noise_colors.push([]);
        noise_greyscales.push([]);

        for(var j = 0; j < canvas_noise.height; j++)
        {
            var greyscale = Math.random() * 255;
            noise_colors[i].push(greyscale);
            noise_greyscales[i].push(greyscale);
            var x = (j*canvas_noise.width + i)*4;
  		    noiseData.data[x]=greyscale;
  		    noiseData.data[x+1]= greyscale;
  		    noiseData.data[x+2]= greyscale;
            noiseData.data[x+3]= 255.0;
        }
    }
    ctx_noise.putImageData(noiseData, 0, 0);
}