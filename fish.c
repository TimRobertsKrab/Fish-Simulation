#include <GL/glut.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// This program will display a simulation of fish motion implementing Couzin's model

#define DEG_TO_RAD 0.017453293
#define PI 3.14159265358979323846

#define MAX_FISH 1000
#define MAX_BOX_EDGE 50
#define MIN_BOX_EDGE 25

struct fish {
    GLfloat position_v[3]; // x y z co-ordinates
    GLfloat direction_v[3]; //x y z co-ordinates unit vector for it's direction
    GLfloat next_direction_v[3]; // the vector direction_v will become.
    int in_ZOR; // identifier for if another fish was in the ZOR
    int in_ZOO; // identifier for if another fish was in the ZOO
    int in_ZOA; // identifier for if another fish was in the ZOA
};

GLfloat  eyex, eyey, eyez;    // Eye point                                     

GLint width = 1280, height = 960;      /* size of window           */
GLfloat box_edge_size = 50.0; // width/2 of box
GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_position0[4];

//Light and colour arrays for objects in scene
GLfloat matSpecular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat matShininess[] = { 50.0 };
GLfloat matSurface[] = { 0.2, 1.0, 0.2, 0.1 };
GLfloat matSurface2[] = { 1.0, 1.0, 1.0, 0.1 };
GLfloat matEmissive[] = { 0.0, 1.0, 0.0, 0.1 };

GLfloat species1_colour[] = { 1.0,0.5,0.0,0.1 }; // Colour of species one (orange)
GLfloat species2_colour[] = { 0.2,0.4,1.0,0.1 }; // Colour of species two (blue)

GLdouble blind_angle = 90.0; // Determines the volume in which a fish can't 'see' other fish within
GLdouble blind_radian_segment; // blind_angle converted to a value that can be used in calculations 
GLdouble turning_angle_spec1 = 5.0; // The turning angle of species one
GLdouble turning_angle_spec2 = 5.0; // The turning angle of species two
GLdouble turning_radian_spec1; // turning_angle_spec1 converted to radians
GLdouble turning_radian_spec2; // turning_angle_spec2 converted to radians

int fish1_count = 100; // Amount of species one fish in the scene
int fish2_count = 100; // Amount of species two fish in the scene

int ZOR_range_spec1[] = { 2,2 }; // Zone of repulsion range for species one {species one, species two}      
int ZOO_range_spec1[] = { 10,0 }; // Zone of orientation range for species one {species one, species two}  
int ZOA_range_spec1[] = { 20,0 }; // Zone of attraction range for species one {species one, species two}  

int ZOR_range_spec2[] = { 2,2 }; // Zone of repulsion range for species two {species one, species two}          
int ZOO_range_spec2[] = { 0,10 }; // Zone of orientation range for species two {species one, species two}
int ZOA_range_spec2[] = { 0,20 }; // Zone of attraction range for species two {species one, species two}

int hard_wall, pause; // Identifier for if the walls "wrap around" and if the simulation if paused
int two_species = 0; // Identifier for if the second species are activated

GLfloat dist_from_scene; // Value used for the camera's viewpoint

struct fish *f1; // Pointer for species one array
struct fish *f2; // Pointer for species two array

                 // Calculates the length of the given vector
GLfloat calculate_magnitude(GLfloat *vector) {
    return(fabs(sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])));
}

// Converts a vector into a unit vector
void normalise_vector(GLfloat *vector) {
    GLfloat m = calculate_magnitude(vector);
    vector[0] /= m;
    vector[1] /= m;
    vector[2] /= m;
}

// Calculates the cross product of 2 vectors
void calculate_cross_prod(GLfloat *v1, GLfloat *v2, GLfloat *cross) {
    cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Calculates the dot product of 2 vectors
GLfloat calculate_dot_prod(GLfloat *v1, GLfloat *v2) {
    int i;
    GLfloat dot_prod = 0.0;
    for (i = 0; i < 3; i++) {
        dot_prod += (GLdouble)(v1[i] * v2[i]);
    }
    return dot_prod;
}

// Calculates the direction vector start_v -> end_v
// Swapping start_v and end_v produces the vector end_v -> start_v)
void calculate_direction_vector(GLfloat *start_v, GLfloat *end_v, GLfloat *dir_v) {
    int i;
    for (i = 0; i < 3; i++)
        dir_v[i] = (end_v[i] - start_v[i]);
}

// Adds the direction vector start_v -> end_v to next_dir_v
void update_direction_vector(GLfloat *start_v, GLfloat *end_v, GLfloat *next_dir_v) {
    int i;
    GLfloat dir_v[3];
    calculate_direction_vector(start_v, end_v, dir_v);
    normalise_vector(dir_v);
    for (i = 0; i < 3; i++) {
        next_dir_v[i] += dir_v[i];
    }
}

// Calculates the angle between the direction vector pos_start -> pos_end and dir_v_A
GLdouble calculate_angle(GLfloat *dir_v_A, GLfloat *dir_v_B) {
    GLdouble ang;
    GLfloat dot_prod = 0.0;

    dot_prod = calculate_dot_prod(dir_v_A, dir_v_B);
    ang = acos(dot_prod / (calculate_magnitude(dir_v_A) * calculate_magnitude(dir_v_B)));
    return ang;
}

// Calculates the distance between 2 position vectors
GLfloat calculate_distance(GLfloat *vec_a, GLfloat *vec_b) {
    GLfloat dir_vec[3];
    int i;
    for (i = 0; i < 3; i++) {
        dir_vec[i] = vec_a[i] - vec_b[i];
    }
    return calculate_magnitude(dir_vec);
}

// Rodrigues' rotation formula
// Rotates a vector around an axis
void rotate_vector(GLfloat *axis_v, GLfloat *dir_v, GLdouble radian) {
    int i;
    GLfloat cross[3];
    calculate_cross_prod(axis_v, dir_v, cross);

    for (i = 0; i < 3; i++) {
        dir_v[i] = (dir_v[i] * cos(radian)) + (cross[i] * sin(radian));
    }
}

//Return random GLfloat within range [-box_edge_size,box_edge_size]
GLfloat generate_box_value() {
    return (((rand() / (double)RAND_MAX)) - 0.5) * box_edge_size * 2.0;
}

// Generates a random vector
void generate_vector(GLfloat *vector) {
    vector[0] = generate_box_value();
    vector[1] = generate_box_value();
    vector[2] = generate_box_value();
}

// Maps a position vector to RGB values    
void calculate_rgb(GLfloat *position, GLfloat *rgb) {
    int i;
    for (i = 0; i < 3; i++) {
        rgb[i] = (position[i] / box_edge_size * 2) + 0.5;
    }
}

// Draws all of the objects in the scene
void draw_scene(void) {
    int x, z, y;

    glEnable(GL_LIGHTING);

    int i;
    for (i = 0; i < fish1_count; i++) {
        if (!two_species) {
            calculate_rgb(f1[i].position_v, matSurface);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matSurface);
        }
        else {
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, species1_colour);
        }
        glPushMatrix();
        glTranslatef(f1[i].position_v[0], f1[i].position_v[1], f1[i].position_v[2]);
        glutSolidSphere(0.5, 20, 20);
        glPopMatrix();
    }
    if (two_species) {
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, species2_colour);
        for (i = 0; i < fish2_count; i++) {
            glPushMatrix();
            glTranslatef(f2[i].position_v[0], f2[i].position_v[1], f2[i].position_v[2]);
            glutSolidSphere(0.5, 20, 20);
            glPopMatrix();
        }
    }

    glMaterialfv(GL_FRONT, GL_DIFFUSE, matSurface2);

    glDepthRange(0.1, 1.0);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_QUADS);
    glNormal3f(0.0, 1.0, 0.0);
    glVertex3f(-box_edge_size, -box_edge_size, box_edge_size);
    glVertex3f(box_edge_size, -box_edge_size, box_edge_size);
    glVertex3f(box_edge_size, -box_edge_size, -box_edge_size);
    glVertex3f(-box_edge_size, -box_edge_size, -box_edge_size);
    glEnd();

    glBegin(GL_QUADS);
    glNormal3f(0.0, 0.0, 1.0);
    glVertex3f(-box_edge_size, -box_edge_size, -box_edge_size);
    glVertex3f(box_edge_size, -box_edge_size, -box_edge_size);
    glVertex3f(box_edge_size, box_edge_size, -box_edge_size);
    glVertex3f(-box_edge_size, box_edge_size, -box_edge_size);
    glEnd();

    glBegin(GL_QUADS);
    glNormal3f(-1.0, 0.0, 0.0);
    glVertex3f(box_edge_size, -box_edge_size, -box_edge_size);
    glVertex3f(box_edge_size, -box_edge_size, box_edge_size);
    glVertex3f(box_edge_size, box_edge_size, box_edge_size);
    glVertex3f(box_edge_size, box_edge_size, -box_edge_size);
    glEnd();

    glBegin(GL_QUADS);
    glNormal3f(0.0, -1.0, 0.0);
    glVertex3f(-box_edge_size, box_edge_size, box_edge_size);
    glVertex3f(-box_edge_size, box_edge_size, -box_edge_size);
    glVertex3f(box_edge_size, box_edge_size, -box_edge_size);
    glVertex3f(box_edge_size, box_edge_size, box_edge_size);
    glEnd();

    glDisable(GL_LIGHTING);
    glDepthRange(0.0, 0.9);
    glColor3f(0.2, 0.2, 0.2);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (x = -box_edge_size; x <= box_edge_size; x += 2) {
        glVertex3f((GLfloat)x, -box_edge_size + 0.01, -box_edge_size);
        glVertex3f((GLfloat)x, -box_edge_size + 0.01, box_edge_size);
    }
    for (z = -box_edge_size; z <= box_edge_size; z += 2) {
        glVertex3f(-box_edge_size, -box_edge_size + 0.01, (GLfloat)z);
        glVertex3f(box_edge_size, -box_edge_size + 0.01, (GLfloat)z);
    }
    glEnd();

    glBegin(GL_LINES);
    for (x = -box_edge_size; x <= box_edge_size; x += 2) {
        glVertex3f((GLfloat)x, -box_edge_size, -box_edge_size + 0.01);
        glVertex3f((GLfloat)x, box_edge_size, -box_edge_size + 0.01);
    }
    for (y = -box_edge_size; y <= box_edge_size; y += 2) {
        glVertex3f(-box_edge_size, (GLfloat)y, -box_edge_size + 0.01);
        glVertex3f(box_edge_size, (GLfloat)y, -box_edge_size + 0.01);
    }
    glEnd();

    glBegin(GL_LINES);
    for (z = -box_edge_size; z <= box_edge_size; z += 2) {
        glVertex3f(box_edge_size - 0.01, -box_edge_size, (GLfloat)z);
        glVertex3f(box_edge_size - 0.01, box_edge_size, (GLfloat)z);
    }
    for (y = -box_edge_size; y <= box_edge_size; y += 2) {
        glVertex3f(box_edge_size - 0.01, y, -box_edge_size);
        glVertex3f(box_edge_size - 0.01, y, box_edge_size);
    }
    glEnd();

    glBegin(GL_LINES);
    for (x = -box_edge_size; x <= box_edge_size; x += 2) {
        glVertex3f((GLfloat)x, box_edge_size - 0.01, -box_edge_size);
        glVertex3f((GLfloat)x, box_edge_size - 0.01, box_edge_size);
    }
    for (z = -box_edge_size; z <= box_edge_size; z += 2) {
        glVertex3f(-box_edge_size, box_edge_size - 0.01, (GLfloat)z);
        glVertex3f(box_edge_size, box_edge_size - 0.01, (GLfloat)z);
    }
    glEnd();
} // draw_scene()


  // Returns the number of digits "number" has
int return_digits(int number) {
    int digits = 0;
    if (number == 0)
        return 1;
    while (number != 0) {
        number /= 10;
        digits++;
    }
    return digits;
}

// Displays the string on the screen
void print_text(char *string, void *font, GLfloat x, GLfloat y) {
    glRasterPos2i(x, y);
    int i;
    for (i = 0; i < strlen(string); i++)
    {
        glutBitmapCharacter(font, string[i]);
    }
}

// Draws text of parameters and input commands
void draw_HUD() {
    glDisable(GL_TEXTURE_2D);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, width, 0.0, height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    void * font = GLUT_BITMAP_9_BY_15;
    char string[20]; // Len = characters plus digits plus '\0'

    GLfloat y_pos = height - 20;

    snprintf(string, 10 + return_digits(fish1_count), "Spec 1#: %d", fish1_count);
    print_text(string, font, 5, y_pos -= 15);

    snprintf(string, 13 + return_digits(turning_angle_spec1), "Turn angle: %f", turning_angle_spec1);
    print_text(string, font, 10, y_pos -= 15);

    snprintf(string, 11 + return_digits(ZOR_range_spec1[0]), "ZOR(1-1): %d", ZOR_range_spec1[0]);
    print_text(string, font, 10, y_pos -= 15);

    snprintf(string, 11 + return_digits(ZOO_range_spec1[0]), "ZOO(1-1): %d", ZOO_range_spec1[0]);
    print_text(string, font, 10, y_pos -= 15);

    snprintf(string, 11 + return_digits(ZOA_range_spec1[0]), "ZOA(1-1): %d", ZOA_range_spec1[0]);
    print_text(string, font, 10, y_pos -= 15);

    if (two_species) {
        snprintf(string, 11 + return_digits(ZOR_range_spec1[1]), "ZOR(1-2): %d", ZOR_range_spec1[1]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOO_range_spec1[1]), "ZOO(1-2): %d", ZOO_range_spec1[1]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOA_range_spec1[1]), "ZOA(1-2): %d", ZOA_range_spec1[1]);
        print_text(string, font, 10, y_pos -= 15);

        y_pos = height - 150;

        snprintf(string, 10 + return_digits(fish2_count), "Spec 2#: %d", fish2_count);
        print_text(string, font, 5, y_pos -= 15);

        snprintf(string, 13 + return_digits(turning_angle_spec2), "Turn angle: %f", turning_angle_spec2);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOR_range_spec2[0]), "ZOR(2-1): %d", ZOR_range_spec2[0]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOO_range_spec2[0]), "ZOO(2-1): %d", ZOO_range_spec2[0]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOA_range_spec2[0]), "ZOA(2-1): %d", ZOA_range_spec2[0]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOR_range_spec2[1]), "ZOR(2-2): %d", ZOR_range_spec2[1]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOO_range_spec2[1]), "ZOO(2-2): %d", ZOO_range_spec2[1]);
        print_text(string, font, 10, y_pos -= 15);

        snprintf(string, 11 + return_digits(ZOA_range_spec2[1]), "ZOA(2-2): %d", ZOA_range_spec2[1]);
        print_text(string, font, 10, y_pos -= 15);
    }
    y_pos = 330;

    print_text("Controls -", font, 5, y_pos -= 15);

    print_text("ZOR(1-1): 'e,r'", font, 10, y_pos -= 15);
    print_text("ZOO(1-1): 'd,f'", font, 10, y_pos -= 15);
    print_text("ZOA(1-1): 'c,v'", font, 10, y_pos -= 15);
    print_text("ZOR(1-2): 'E,R'", font, 10, y_pos -= 15);
    print_text("ZOO(1-2): 'D,F'", font, 10, y_pos -= 15);
    print_text("ZOA(1-2): 'C,V'", font, 10, y_pos -= 15);
    print_text("Turn angle 1: ', .'", font, 10, y_pos -= 15);

    print_text("ZOR(2-1): 'y,u'", font, 10, y_pos -= 15);
    print_text("ZOO(2-1): 'h,j'", font, 10, y_pos -= 15);
    print_text("ZOA(2-1): 'n,m'", font, 10, y_pos -= 15);
    print_text("ZOR(2-2): 'Y,U'", font, 10, y_pos -= 15);
    print_text("ZOO(2-2): 'H,J'", font, 10, y_pos -= 15);
    print_text("ZOA(2-2): 'N,M'", font, 10, y_pos -= 15);
    print_text("Turn angle 2: '< >'", font, 10, y_pos -= 15);
    print_text("Change view: Left & Right", font, 10, y_pos -= 15);
    print_text("Change size: Up & Down", font, 10, y_pos -= 15);
    print_text("Restart: 'q'", font, 10, y_pos -= 15);
    print_text("Toggle walls: 'a'", font, 10, y_pos -= 15);
    print_text("Toggle species: 'z'", font, 10, y_pos -= 15);
    print_text("Pause: 'p'", font, 10, y_pos -= 15);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_TEXTURE_2D);
}

// Manages material properties and drawing on the screen
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
    glMaterialfv(GL_FRONT, GL_SHININESS, matShininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matSurface);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matSurface2);

    glLoadIdentity();
    draw_HUD();
    gluLookAt(eyex, eyey, eyez, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    draw_scene();
    glutSwapBuffers();
}

// Sets all the dimensions of a vector to zero
void initialise_vector(GLfloat *vector) {
    vector[0] = vector[1] = vector[2] = 0;
}

// Checks if a vector is a zero vector
int is_zero_vector(GLfloat *vector) {
    if (vector[0] == 0 && vector[1] == 0 && vector[2] == 0)
        return 1;

    return 0;
}

// Rotates the direction vector closer to it's next direction vector, alters the position vector
// and manages wall collision
void move_fish(struct fish *f, int fish_count, GLdouble radian) {
    int i, j;
    GLfloat normal_v[3], safety_v[3];


    for (i = 0; i < fish_count; i++) {
        if (!(is_zero_vector(f[i].next_direction_v))) {
            if (calculate_angle(f[i].direction_v, f[i].next_direction_v) <= radian) {
                for (j = 0; j < 3; j++) {
                    f[i].direction_v[j] = f[i].next_direction_v[j];
                }
            }
            else {
                calculate_cross_prod(f[i].direction_v, f[i].next_direction_v, normal_v);
                // If direction_v == -(next_direction_v) the normal will be the zero vector
                // This causes a chain reaction in which the position will be set to NaN.
                if (!(is_zero_vector(normal_v))) {
                    normalise_vector(normal_v);
                    rotate_vector(normal_v, f[i].direction_v, radian);
                }
                else {
                    do {
                        generate_vector(safety_v);
                        normalise_vector(safety_v);
                    } while (f[i].direction_v[0] == safety_v[0] && f[i].direction_v[1] == safety_v[1] && f[i].direction_v[1] == safety_v[1]);
                    calculate_cross_prod(f[i].direction_v, safety_v, normal_v);
                    rotate_vector(normal_v, f[i].direction_v, radian);
                }
            }
        }
        for (j = 0; j < 3; j++) {
            f[i].position_v[j] += f[i].direction_v[j];
            if (fabs(f[i].position_v[j]) > box_edge_size) {
                f[i].position_v[j] *= (box_edge_size / fabs(f[i].position_v[j]));
                if (hard_wall)
                    f[i].direction_v[j] *= -1.0;
                else
                    f[i].position_v[j] *= -1.0;
            }
        }

    }
}

// Determines the next direction vector with regards to the zone of repulsion
void update_in_ZOR(struct fish *fish1, struct fish *fish2, int fish_count, int zor) {
    int i;
    GLfloat vector[3];

    for (i = 0; i < fish_count; i++) {
        if (fish1 != &fish2[i]) {
            calculate_direction_vector(fish1->position_v, fish2[i].position_v, vector);
            if (calculate_distance(fish1->position_v, fish2[i].position_v) < zor &&
                calculate_angle(fish1->direction_v, vector) < blind_radian_segment) {
                fish1->in_ZOR = 1;
                update_direction_vector(fish2[i].position_v, fish1->position_v, fish1->next_direction_v);
            }
        }
    }
}

// Determines the next direction vector with regards to the zone of orientation 
// and zone of attraction
void update_in_ZOO_ZOA(struct fish *fish1, struct fish *fish2, int fish_count, int zoo, int zoa) {
    int i, j;
    GLfloat dist;
    GLfloat vector[3];

    for (i = 0; i < fish_count; i++) {
        if (fish1 != &fish2[i]) {
            calculate_direction_vector(fish1->position_v, fish2[i].position_v, vector);
            if (calculate_angle(fish1->direction_v, vector) < blind_radian_segment) {
                dist = calculate_distance(fish1->position_v, fish2[i].position_v);
                if (dist < zoo) {
                    fish1->in_ZOO = 1;
                    for (j = 0; j < 3; j++) {
                        fish1->next_direction_v[j] += fish2[i].direction_v[j];
                    }
                }
                else if (dist >= zoo && dist < zoa) {
                    fish1->in_ZOA = 1;
                    update_direction_vector(fish1->position_v, fish2[i].position_v, fish1->next_direction_v);
                }
            }
        }
    }
}


// Updates the positions and directions of the fish.
void update_fish(void) {
    int i, j;

    if (!pause) {
        // Alter fish positions
        move_fish(f1, fish1_count, turning_radian_spec1);
        if (two_species) {
            move_fish(f2, fish2_count, turning_radian_spec2);
        }
        // Alter next_direction vectors of species 1.
        for (i = 0; i < fish1_count; i++) {
            initialise_vector(f1[i].next_direction_v);
            update_in_ZOR(&f1[i], f1, fish1_count, ZOR_range_spec1[0]);
            if (two_species) {
                update_in_ZOR(&f1[i], f2, fish2_count, ZOR_range_spec1[1]);
            }
            // Only do ZOO,ZOA work if no fish were in the ZOR
            if (!(f1[i].in_ZOR)) {
                update_in_ZOO_ZOA(&f1[i], f1, fish1_count, ZOO_range_spec1[0], ZOA_range_spec1[0]);
                if (two_species) {
                    update_in_ZOO_ZOA(&f1[i], f2, fish2_count, ZOO_range_spec1[1], ZOA_range_spec1[1]);
                }
            }
        }

        // Alter next_direction vectors of species 2.
        if (two_species) {
            for (i = 0; i < fish2_count; i++) {
                initialise_vector(f2[i].next_direction_v);
                update_in_ZOR(&f2[i], f1, fish1_count, ZOR_range_spec2[0]);
                update_in_ZOR(&f2[i], f2, fish2_count, ZOR_range_spec2[1]);
                // Only do ZOO,ZOA work if no fish were in the ZOR
                if (!(f2[i].in_ZOR)) {
                    update_in_ZOO_ZOA(&f2[i], f1, fish1_count, ZOO_range_spec2[0], ZOA_range_spec2[0]);
                    update_in_ZOO_ZOA(&f2[i], f2, fish2_count, ZOO_range_spec2[1], ZOA_range_spec2[1]);
                }
            }
        }

        for (i = 0; i < fish1_count; i++) {
            if (f1[i].in_ZOO) {
                for (j = 0; j < 3; j++) {
                    f1[i].next_direction_v[j] += f1[i].direction_v[j];
                }
            }
            if (f1[i].in_ZOR || f1[i].in_ZOO || f1[i].in_ZOA) {
                normalise_vector(f1[i].next_direction_v);
                f1[i].in_ZOR = 0;
                f1[i].in_ZOO = 0;
                f1[i].in_ZOA = 0;
            }
        }
        if (two_species) {
            for (i = 0; i < fish2_count; i++) {
                if (f2[i].in_ZOO) {
                    for (j = 0; j < 3; j++) {
                        f2[i].next_direction_v[j] += f2[i].direction_v[j];
                    }
                }
                if (f2[i].in_ZOR || f2[i].in_ZOO || f2[i].in_ZOA) {
                    normalise_vector(f2[i].next_direction_v);
                    f2[i].in_ZOR = 0;
                    f2[i].in_ZOO = 0;
                    f2[i].in_ZOA = 0;
                }
            }
        }
    }
    glutPostRedisplay();
}

// Manages the window when it is reshaped    
void reshape(int w, int h) {
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50, (GLfloat)w / (GLfloat)h, 0.1, 1000);
    glMatrixMode(GL_MODELVIEW);
    width = w;
    height = h;
}

// Changes the viewpoint and size of the "box"    
void cursor_keys(int key, int x, int y) {
    switch (key) {
    case GLUT_KEY_RIGHT:
        if (eyez == 0.0) {
            eyez = box_edge_size + 65.0;
            eyex = -box_edge_size - 65.0;
        }
        else {
            eyex = 0.0;
            eyez = dist_from_scene;
        }
        break;
    case GLUT_KEY_LEFT:
        if (eyex == 0.0) {
            eyez = box_edge_size + 65.0;
            eyex = -box_edge_size - 65.0;
        }
        else {
            eyez = 0.0;
            eyex = dist_from_scene * -1;
        }
        break;
    case GLUT_KEY_UP:
        if (box_edge_size < MAX_BOX_EDGE)
            box_edge_size++;
        break;
    case GLUT_KEY_DOWN:
        if (box_edge_size > MIN_BOX_EDGE)
            box_edge_size--;
        break;
    }
}

// Allocates memory and initialises fish variables.
void init(void) {
    light_position0[0] = -box_edge_size;
    light_position0[1] = light_position0[3] = 0.0;
    light_position0[2] = box_edge_size;
    glClearColor(1.0, 1.0, 1.0, 0.0);   /* Define background colour */
    int i;
    // Allocate memory for fish
    while (f1 == NULL)
        f1 = (struct fish*)malloc(sizeof(struct fish) * MAX_FISH);
    while (f2 == NULL)
        f2 = (struct fish*)malloc(sizeof(struct fish) * MAX_FISH);
    // Initialise the fish
    for (i = 0; i < MAX_FISH; i++) {
        generate_vector(f1[i].position_v);
        generate_vector(f1[i].direction_v);
        normalise_vector(f1[i].direction_v);
        initialise_vector(f1[i].next_direction_v);
        f1[i].in_ZOR = 0;
        f1[i].in_ZOO = 0;
        f1[i].in_ZOA = 0;
        generate_vector(f2[i].position_v);
        generate_vector(f2[i].direction_v);
        normalise_vector(f2[i].direction_v);
        initialise_vector(f2[i].next_direction_v);
        f2[i].in_ZOR = 0;
        f2[i].in_ZOO = 0;
        f2[i].in_ZOA = 0;
    }
    blind_radian_segment = PI - (blind_angle * DEG_TO_RAD * 0.5);
    hard_wall = 1;
    pause = 0;
    eyex = -box_edge_size - 65.0;
    eyey = 0.0;
    eyez = box_edge_size + 65.0;
    dist_from_scene = sqrt(fabs(eyex) * fabs(eyex) + fabs(eyez) * fabs(eyez));

    // Set up lighting 
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
    glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
}

// Allows zone ranges, fish counts, wall state and species state to be altered.
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 27:
        free(f1);
        free(f2);
        exit(0);
        break;
    case 'q':
        init();
        break;
    case 'a':
        hard_wall = !hard_wall;
        break;
    case 'z':
        two_species = !two_species;
        break;
    case '[':
        if (fish1_count > 0) {
            fish1_count--;
        }
        break;
    case ']':
        if (fish1_count < MAX_FISH)
            fish1_count++;
        break;
    case '{':
        if (two_species)
            if (fish2_count > 0)
                fish2_count--;
        break;
    case '}':
        if (two_species)
            if (fish2_count < MAX_FISH)
                fish2_count++;
        break;
    case 'e':
        if (ZOR_range_spec1[0] > 0)
            ZOR_range_spec1[0]--;
        break;
    case 'r':
        if (ZOR_range_spec1[0] < 2 * box_edge_size)
            ZOR_range_spec1[0]++;
        break;
    case 'd':
        if (ZOO_range_spec1[0] > 0)
            ZOO_range_spec1[0]--;
        break;
    case 'f':
        if (ZOO_range_spec1[0] < 2 * box_edge_size)
            ZOO_range_spec1[0]++;
        break;
    case 'c':
        if (ZOA_range_spec1[0] > 0)
            ZOA_range_spec1[0]--;
        break;
    case 'v':
        if (ZOA_range_spec1[0] < 2 * box_edge_size)
            ZOA_range_spec1[0]++;
        break;
    case 'E':
        if (ZOR_range_spec1[1] > 0)
            ZOR_range_spec1[1]--;
        break;
    case 'R':
        if (ZOR_range_spec1[1] < 2 * box_edge_size)
            ZOR_range_spec1[1]++;
        break;
    case 'D':
        if (ZOO_range_spec1[1] > 0)
            ZOO_range_spec1[1]--;
        break;
    case 'F':
        if (ZOO_range_spec1[1] <  2 * box_edge_size)
            ZOO_range_spec1[1]++;
        break;
    case 'C':
        if (ZOA_range_spec1[1] > 0)
            ZOA_range_spec1[1]--;
        break;
    case 'V':
        if (ZOA_range_spec1[1] < 2 * box_edge_size)
            ZOA_range_spec1[1]++;
        break;
    case ',':
        if (turning_angle_spec1 > 1)
            turning_radian_spec1 = (turning_angle_spec1--) * DEG_TO_RAD;
        break;
    case '.':
        if (turning_angle_spec1 < 10)
            turning_radian_spec1 = (turning_angle_spec1++) * DEG_TO_RAD;
        break;
    case 'y':
        if (ZOR_range_spec2[0] > 0)
            ZOR_range_spec2[0]--;
        break;
    case 'u':
        if (ZOR_range_spec2[0] < 2 * box_edge_size)
            ZOR_range_spec2[0]++;
        break;
    case 'h':
        if (ZOO_range_spec2[0] > 0)
            ZOO_range_spec2[0]--;
        break;
    case  'j':
        if (ZOO_range_spec2[0] < 2 * box_edge_size)
            ZOO_range_spec2[0]++;
        break;
    case 'n':
        if (ZOA_range_spec2[0] > 0)
            ZOA_range_spec2[0]--;
        break;
    case 'm':
        if (ZOA_range_spec2[0] < 2 * box_edge_size)
            ZOA_range_spec2[0]++;
        break;
    case 'Y':
        if (ZOR_range_spec2[1] > 0)
            ZOR_range_spec2[1]--;
        break;
    case 'U':
        if (ZOR_range_spec2[1] < 2 * box_edge_size)
            ZOR_range_spec2[1]++;
        break;
    case 'H':
        if (ZOO_range_spec2[1] > 0)
            ZOO_range_spec2[1]--;
        break;
    case 'J':
        if (ZOO_range_spec2[1] <  2 * box_edge_size)
            ZOO_range_spec2[1]++;
        break;
    case 'N':
        if (ZOA_range_spec2[1] > 0)
            ZOA_range_spec2[1]--;
        break;
    case 'M':
        if (ZOA_range_spec2[1] < 2 * box_edge_size)
            ZOA_range_spec2[1]++;
        break;
    case '<':
        if (turning_angle_spec2 > 1)
            turning_radian_spec2 = (turning_angle_spec2--) * DEG_TO_RAD;
        break;
    case '>':
        if (turning_angle_spec2 < 10)
            turning_radian_spec2 = (turning_angle_spec2++) * DEG_TO_RAD;
        break;
    case 'p':
        pause = !pause;
    }
}

// Main method    
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("Simulation of fish motion");
   // glutFullScreen();
    turning_radian_spec1 = turning_angle_spec1 * DEG_TO_RAD;
    turning_radian_spec2 = turning_angle_spec2 * DEG_TO_RAD;
    init();
    glutDisplayFunc(display);
    glutIdleFunc(update_fish);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(cursor_keys);
    glutMainLoop();
    return 0;
}
