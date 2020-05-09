/*
credits: Matthias Hasler @drhasler
He is very good at programming and he is nice
also he made a thing to eat cats ...ew
*/
#include<bits/stdc++.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
mesh eat_da_mesh(){
    std::ifstream is("meshy_bois/cat.obj");
    if (!is.is_open()) { std::cerr << "couldnt open it\n"; exit(1); }
    std::string s; mesh m;
    while (is >> s) {
        if (s[0] == 'v') {
            double x,y,z;
            is >> x >> y >> z;
            ( s[1]=='n' ? m.norm : s[1]=='t' ? m.uv : m.vertices )
                .push_back(Vector{x,y,z});
        } else if (s[0] == 'f') {
            char c; triangle_indices t;
            for (int i=0; i<3; i++)
                is >> t.vtx[i] >> c >> t.uv[i] >> c >> t.norm[i],
            t.vtx[i] --,
            t.uv[i]  --, //some autists decided to 1-index this stuff so we counter autism with autism :)
            t.norm[i]--;
            m.triangs.push_back(t);
        }
        is.ignore(256,'\n');
    }
    int nc;
    m.image = stbi_load("meshy_bois/cat_diff.png",
            &m.width, &m.height, &nc, 3);
    //m.shift(Vector{0,-20,0});
    //m.scale(1.2);
    //m.rot(-2.5,0,0);
    //m.make_bvh();
    return m;
}
