#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

float pi = 3.14159265359;
float degToRad = 3.14159265359 / 180.0;

vec2 smin( vec2 a, vec2 b, float k )
{
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    if(h < 0.5) {
        return vec2(mix( b.x, a.x, h ) - k*h*(1.0-h), b.y);
    } else {
        return vec2(mix( b.x, a.x, h ) - k*h*(1.0-h), a.y);
    }
}

vec3 elongate( vec3 p, vec3 h )
{
    vec3 q = p - clamp( p, -h, h );
    return q;
}
float cappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;
    float x = sqrt( papa - paba*paba*baba );
    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;
    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
    float cbx = x-ra - f*rba;
    float cby = paba - f;
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

float roundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;
    
    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot( pa*l2 - ba*y, pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

vec2 smax( vec2 a, vec2 b, float k )
{
    float h = max(k-abs(a.x-b.x),0.0);
    return (a.x > b.x ? vec2((a.x + h*h*0.25/k), a.y) : vec2((b.x + h*h*0.25/k), b.y));
}


// Ra: radius rb: roundedness h: height
float roundedCylinder( vec3 p, float ra, float rb, float h )
{
  vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

float cappedTorus(vec3 p, vec2 sc, float ra, float rb)
{
  p.x = abs(p.x);
  float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
  return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}

float box( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}


float sphere(vec3 p, float s)
{
    return length(p) - s;
}

float ellipsoid( vec3 p, vec3 r )
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float stick(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    vec3 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return  length( pa - ba*h ) - mix(r1,r2,h*h*(3.0-2.0*h));
}


vec3 bend( vec3 p, float k )
{
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xy,p.z);
    return q;
}


vec3 twist( vec3 p, float k)
{
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xz,p.y);
    return q;
}

float onion(float sdf, float thickness)
{
    return abs(sdf) - thickness;
}

float roundBox( vec3 p, vec3 b, float r )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

// https://iquilezles.org/www/articles/noacos/noacos.htm
mat3 rotationAxisAngle( vec3 v, float a )
{
    float si = sin( a );
    float co = cos( a );
    float ic = 1.0f - co;

    return mat3( v.x*v.x*ic + co,       v.y*v.x*ic - si*v.z,    v.z*v.x*ic + si*v.y,
                   v.x*v.y*ic + si*v.z,   v.y*v.y*ic + co,        v.z*v.y*ic - si*v.x,
                   v.x*v.z*ic - si*v.y,   v.y*v.z*ic + si*v.x,    v.z*v.z*ic + co );
}

float plane( vec3 p, vec3 n, float h )
{
  // n must be normalized
  return dot(p,n) + h;
}

float getBias(float time, float bias)
{
    return (time / ((((1.0/bias) - 2.0)*(1.0 - time))+1.0));
}

float getGain(float time,float gain)
{
  if(time < 0.5)
    return getBias(time * 2.0,gain)/2.0;
  else
    return getBias(time * 2.0 - 1.0,1.0 - gain)/2.0 + 0.5;
}

vec2 map(vec3 p)
{
//    if (abs(p.z) > 60.0 || abs(p.x) > 30.0 || abs(p.y) > 30.0)
//    {
//        return vec2(10000.0, 0.0);
//    }
    vec2 face = vec2(box(p - vec3(0.0,-7.4,0.0), vec3(9.0,14.0,9.0)), 0.1);
    
    if(face.x < 0.0001) {
    vec3 symP = vec3(abs(p.x), p.yz);
    //TODO: calculate these matrices on CPU or vert shader
    mat3 rot = rotationAxisAngle(normalize(vec3(1.0, 0.0, 0.0)), degToRad * 10.0);
    face = vec2(sphere(p + vec3(0.0, -0.8, 0.0), 1.2), 0.0);
    vec3 rotP = rot * (p - vec3(0.0,0.9,-0.63));
    face = smin(face, vec2(ellipsoid(p - vec3(0.0,0.0,-0.4), vec3(1.15, 1.3, 1.0)), 0.0), 0.2);
    face = smin(face, vec2(roundBox(rotP, vec3(0.6, 0.2, 0.55), 0.22), 0.0), 0.6);
    // Chin
    face = smin(face, vec2(sphere(p + vec3(0.0, 0.98, 0.89), 0.33), 0.0), 0.3);
    // Eye socket
    face = smax(face, -vec2(sphere(symP + vec3(-0.5, -0.46, 1.2), 0.15), 0.0), 0.4);

    vec3 eyePos = vec3(-0.05, 0.03, 0.15);
    //eyelids
    //top eyelid
    vec3 blinkOffset = vec3(0.0, 0.1 + -0.2 * abs(sin(0.05 * u_Time)), 0.0);
    //blinkOffset.y = 0.1;
    face = smin(face, vec2(stick(symP, vec3(0.4, 0.52, -1.3) + eyePos, vec3(0.75, 0.53, -1.28) + eyePos, 0.075, 0.075), 0.0), 0.07);
    
    face = smin(face, vec2(stick(symP, vec3(0.4, 0.52, -1.3) + eyePos + blinkOffset, vec3(0.7, 0.53, -1.28) + eyePos + blinkOffset, 0.075, 0.075), 0.0), 0.05);

    // low eyelid
    face = smin(face, vec2(stick(symP, vec3(0.76, 0.31, -1.15) + eyePos, vec3(0.3, 0.27, -1.37) + eyePos, 0.05, 0.05), 0.0), 0.08);
    face = smin(face, vec2(sphere(symP + vec3(-0.45, -0.04, 0.88), 0.36), 0.0), 0.15);

    // Eyeball
    vec3 eyeballPos = symP + vec3(-0.5, -0.44, 0.99);
    vec2 eyeball = vec2(sphere(eyeballPos, 0.19), 3.0);
    float pupilDist = length(eyeballPos.xy);
    if(pupilDist < 0.08 && pupilDist > 0.04) {
        eyeball.y = 5.0;
    }
    face = smin(face, eyeball, 0.0);
    
    // Cheekbones
    face = smin(face, vec2(stick(symP, vec3(0.92, 0.22, -0.64), vec3(0.52, -0.38, -0.84), 0.2, 0.25), 0.0), 0.24);
    // Hair Loop
    vec3 bentSymP = bend(symP, -0.7);
    rot = rotationAxisAngle(normalize(vec3(0.8, 0.0, 1.0)), -degToRad * 60.0);
    vec3 rotSymP = rot * (symP + vec3(-1.7, -2.4, -0.3));
    rotSymP = bend(rotSymP, 0.05);
    rotSymP += 0.05 * sin(length(rotSymP * 4.0)) + 0.05 * cos(p * 3.5);
    vec2 hair = vec2(roundedCylinder(rotSymP, 1.1, 0.3, 0.1), 1.0);
    
    //scalp hair
    hair = smin(hair, vec2(sphere(symP + vec3(-0.3, -1.2, -0.3), 1.33), 1.0), 0.4);
    face = smin(face, hair, 0.1);
    vec3 twistPos = p;
    twistPos.xz += 0.03 * abs(sin(30.0 * p.y));
    vec2 neck = vec2(stick(twistPos, vec3(0.0, -0.4, 0.4), vec3(0.0, -2.8, 0.8), 0.8, 1.0), 2.0);
    face = smin(face, neck, 0.1);
    
    //lips
    //vec3 bent = bend(symP + vec3(-0.1, 0.5, 1.14), -5.8);
    vec3 lipPos = p + vec3(0.0, 0.49, 1.18);
    float xStretch = lipPos.x * 9.0;
    
    // Quadratic curve for lips
    float offset = -0.03 -0.08 * ((xStretch * xStretch * xStretch * xStretch) - 1.7 * (xStretch * xStretch));
    vec2 lip = vec2(roundBox(lipPos + vec3(0.0, -offset * 0.5, 0.0), vec3(0.16, 0.5 * offset + 0.01  * clamp(1.0 - (abs(lipPos.x * 200.0 * lipPos.x)), 0.0, 1.0), 0.1), 0.1), 0.0);
    
    //lipPos.y -= 0.01;
    lipPos = p + vec3(0.0, 0.52, 1.19);

    vec2 upperLip = vec2(roundBox(lipPos, vec3(0.2, 0.02, max(0.1 - abs(0.12 * p.x), 0.02)), 0.05), 0.0);
    lip = smin(lip, upperLip, 0.02);
    
    lipPos = p + vec3(0.0, 0.64, 1.14);
    float lipFunc = -10.0 * lipPos.x * lipPos.x + 0.4;
    if(lipFunc > lipPos.y) {
        lip.y = 4.0;
    }

    lipPos = p + vec3(0.0, 0.64, 1.14);
    lipFunc = 10.0 * lipPos.x * lipPos.x - 0.2;
    
    lipPos.y += 0.03 * (cos(13.0 * lipPos.x));
    
    vec2 lowerLip = vec2(roundBox(lipPos, vec3(0.13, 0.003, 0.1), 0.1), 0.0);
    
    if(lipFunc < lipPos.y) {
        lowerLip.y = 4.0;
    }
    
    lip = smin (lip, lowerLip, 0.002);
    face = smin(face, lip, 0.01);

    vec3 nosePos = vec3(0.0, -0.1, 0.0);
    vec2 nose = vec2(sphere(p + vec3(0.0,-0.02,1.55) - nosePos, 0.135), 0.0);
    
    // Nose bridge
    nose = smin(nose, vec2(stick(p, vec3(0.0, 0.5, -1.3) + nosePos, vec3(0.0, 0.14, -1.42) + nosePos, 0.11, 0.13), 0.0), 0.12);
    
    // Nose sides
    nose = smin(nose, vec2(stick(symP, vec3(0.1, 0.0, -1.55) + nosePos, vec3(0.13, -0.1, -1.45) + nosePos, 0.06, 0.06), 0.0), 0.09);
    face = smin(face, nose, 0.07);

    //body
    vec3 robeOffset = vec3(0.0, -0.4, -0.6);
    vec3 robeP = (p - vec3(0.0, -3.0, 0.0) + robeOffset);
    float fold = 0.1 * (getGain(0.5 * (sin(0.1 * u_Time + p.x) + 1.0), 0.4));
    robeP.z += fold;
    robeP = elongate(robeP, vec3(1.0, 0.0, 0.0));

    vec2 robe = vec2(cappedCone(robeP, vec3(0.0, 2.6, 0.4), vec3(0.0, -3.4, -0.5), 2.4, 1.4), 5.0);
    if (robe.x < 0.0) {
        robe.y = 4.0;
    }
    
    vec3 foldP = p - vec3(0.0, -13.4, 0.9);
    foldP.z += 0.2 * (getGain(0.5 * (sin(0.1 * u_Time + p.x) + 1.0), 0.4) - 0.5) * 2.0;

    vec2 bottomRobe = vec2(roundBox(foldP, vec3(1.4, 6.9, 0.7), 1.6), 5.0);
    
    bottomRobe.x += 0.001 * (0.5 * sin(p.x * 2.0) + sin(p.y * 1.0) + 2.0 * sin(p.z * 0.5));
    bottomRobe.x = onion(bottomRobe.x, abs(bottomRobe.x) * 0.1);
    bottomRobe = smin(bottomRobe, vec2(cappedCone(p, vec3(0.0, -16.6, 0.9), vec3(0.0, -22.9, -0.5), 2.1, 4.8), 5.0), 0.9);

    foldP = elongate(p, vec3(3.0, 0.0, 0.0));
    foldP += 0.3 * (sin(0.1 * u_Time + p.x) + 1.0);
    //foldP.y += 0.3 * sin(0.1 * u_Time + p.x);
    vec2 robeFolds = vec2(cappedCone(foldP, vec3(0.0, -18.9, 0.9), vec3(0.0, -22.9, -0.5), 0.2, 6.8), 5.0);
    bottomRobe = smin(bottomRobe, robeFolds, 0.3);

    robe = smin(robe, bottomRobe, 0.1);
    vec2 insideRobe = robe.xy;

    //body
    vec2 body = smax(insideRobe, vec2(ellipsoid(p - vec3(0.0,-3.7,0.9), vec3(4.0, 1.4, 1.9)), 0.0), 0.01);
    face = smin(face, body, 0.1);
    robe.x = onion(robe.x, 0.03);
    
    vec2 belt = vec2(onion(robe.x, 0.1), 0.0);
    belt = smax(belt, vec2(roundBox(p - vec3(0.0, -4.5, 0.0), vec3(2.5, 0.9, 3.0), 0.01), 0.0), 0.01);
    robe = smin(robe, belt, 0.001);
    
    vec3 longP = p - vec3(4.0, -6.0, -0.7);//elongate(p - vec3(4.0, -6.0, -0.7), vec3(0.0, 0.4, 0.0));
    vec2 bicep = vec2(stick(symP, vec3(2.3, -3.0, 0.4), vec3(4.0, -5.6, -1.3), 0.9, 1.0), 5.0);
    bicep = smax(bicep, -insideRobe, 0.2);
    //robe = smin(robe, bicep, 0.2);
    vec3 sleeveP = symP - vec3(4.0, -6.0, -0.7);
    vec2 sleeve = vec2(cappedCone(sleeveP, vec3(0.0, 0.0, -0.6), vec3(-4.4, -0.8, -4.0), 0.9, 1.8), 5.0);
    if (sleeve.x < 0.0) {
        sleeve.y = 4.0;
    }
    
    sleeve.x = onion(sleeve.x, 0.1);
    rot = rotationAxisAngle(vec3(0.0, 1.0, 0.0), -degToRad * 10.0);
    
    // Cut out sleeves using box
    rotP = symP - vec3(0.0, -4.0, -4.5);
    rotP.x += 0.05 * getBias((sin(0.3 * u_Time + 4.0 * rotP.y) + 1.0) * 0.5, 0.6);
    
    rotP = rot * (rotP);
    sleeve = smax(sleeve, vec2(-roundBox(rotP, vec3(0.2, 5.0, 4.0), 0.1), 3.0), 0.1);
        
    longP = p - vec3(-4.0, -6.0, -0.7);
    
    vec2 lWrist = vec2(cappedCone(p - vec3(4.0, -5.2, -1.5), vec3(-3.0, -0.4, -1.7), vec3(-5.4, -0.4, -2.5), 0.3, 0.4), 0.0);

    vec2 rWrist = vec2(cappedCone(p - vec3(-4.0, -5.3, -1.6), vec3(3.0, -0.6, -1.8), vec3(5.4, -0.4, -1.5), 0.4, 0.4), 0.0);

    sleeve = smin(lWrist, sleeve, 0.01);
    sleeve = smin(rWrist, sleeve, 0.01);

    vec2 arm = smin(bicep, sleeve, 0.5);
    robe = smin(robe, arm, 0.01);
    // Cut off top of robe
    vec3 clothP = p - vec3(0.0, 2.45, 0.0);
    clothP.y += 0.03 * sin(0.3 * u_Time + 4.0 * clothP.x);
    robe = smax(robe, vec2(-roundBox(clothP, vec3(4.0, 2.9, 4.5), 0.01), 0.0), 0.01);
    rot = rotationAxisAngle(vec3(0.0, 0.0, 1.0), degToRad *  45.0);
    rotP = rot * (p - vec3(0.0, 1.0, -1.2) + robeOffset);
    //Cut out vneck
    robe = smax(robe, vec2(-roundBox(rotP, vec3(3.5, 3.5, 1.0), 0.01), 0.0), 0.01);
    face = smin(face, robe, 0.01);
    }
    return face;
}

vec3 calcNormals(vec3 p)
{
    float epsilon = 0.00001;
    return normalize(vec3(map(p + vec3(epsilon, 0.0, 0.0)).x - map(p - vec3(epsilon, 0.0, 0.0)).x, map(p + vec3(0.0, epsilon, 0.0)).x - map(p - vec3(0.0, epsilon, 0.0)).x, map(p + vec3(0.0, 0.0, epsilon)).x - map(p - vec3(0.0, 0.0, epsilon)).x));
    
}

vec4 raycast(vec3 origin, vec3 dir, int maxSteps)
{
    float t = 0.0;

    for(int i = 0; i < maxSteps; ++i)
    {
        vec3 p = origin + t * dir;
        vec2 dist = map(p);

        if (abs(dist.x) < 0.0001) {
            return vec4(p, dist.y);
        }
        
        t += dist.x;
        if(t > 40.0) {
            return vec4(0.0, 0.0, 0.0, -100.0);

        }
    }
    
    return vec4(0.0, 0.0, 0.0, -100.0);
}

void main() {
    float fov = 22.5f;
    float len = distance(u_Ref, u_Eye);
    vec3 look = normalize(u_Ref - u_Eye);
    vec3 right = normalize(cross(look, u_Up));
    float aspect = u_Dimensions.x / u_Dimensions.y;
    vec3 v = u_Up * len * tan(fov);
    vec3 h = right * len * aspect * tan(fov);

    vec3 p = u_Ref + fs_Pos.x * h + fs_Pos.y * v;
    vec3 dir = normalize(p - u_Eye);
    
    vec3 lightPos = vec3(3.0, 10.5, -20.0);
    
    vec4 isect = raycast(u_Eye, dir, 128);
    float a = getBias((fs_Pos.y + 1.0) * 0.5, 0.3);
    vec3 col = mix(vec3(147.0,157.0,182.0) / 255.0, vec3(226.0,234.0,236.0) / 255.0, a);

    if(isect.w >= 0.0)
    {
        vec3 lightVec = normalize(lightPos - isect.xyz);
        vec3 normal = calcNormals(isect.xyz);
        float ambient = 0.2;
        float diffuse = clamp(dot(normal, lightVec), 0.0, 1.0) + ambient;
        
        vec4 shadow = raycast(isect.xyz + normal * 0.001, lightVec, 120);
        //vec4 shadow = vec4(0.0);
        if (shadow.w >= 0.0) {
            diffuse = ambient;
        }
        
        if(isect.w == 0.0)
        {
            col = diffuse * vec3(0.8, 0.7, 0.76);
        } else if (isect.w <= 1.0) {
            col = diffuse * vec3(0.2, 0.2, 0.2);

        } else if (isect.w <= 2.0) {
            col = diffuse * vec3(0.8, 0.6, 0.1);
        } else if (isect.w <= 3.0) {
            col = diffuse * vec3(0.2, 0.2, 0.2);

        } else if (isect.w <= 4.0) {
            col = diffuse * vec3(0.8, 0.2, 0.4);

        } else if (isect.w <= 5.0) {
            col = diffuse * vec3(177.0,180.0,210.0) / 255.0;

        }


    }
        
    out_Col = vec4(col, 1.0);
}
