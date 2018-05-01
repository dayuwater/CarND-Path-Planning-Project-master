#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{
    
    double closestLen = 100000; //large number
    int closestWaypoint = 0;
    
    for(int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x,y,map_x,map_y);
        if(dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }
        
    }
    
    return closestWaypoint;
    
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    
    int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
    
    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];
    
    double heading = atan2((map_y-y),(map_x-x));
    
    double angle = fabs(theta-heading);
    angle = min(2*pi() - angle, angle);
    
    if(angle > pi()/4)
    {
        closestWaypoint++;
        if (closestWaypoint == maps_x.size())
        {
            closestWaypoint = 0;
        }
    }
    
    return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);
    
    int prev_wp;
    prev_wp = next_wp-1;
    if(next_wp == 0)
    {
        prev_wp  = maps_x.size()-1;
    }
    
    double n_x = maps_x[next_wp]-maps_x[prev_wp];
    double n_y = maps_y[next_wp]-maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];
    
    // find the projection of x onto n
    double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
    double proj_x = proj_norm*n_x;
    double proj_y = proj_norm*n_y;
    
    double frenet_d = distance(x_x,x_y,proj_x,proj_y);
    
    //see if d value is positive or negative by comparing it to a center point
    
    double center_x = 1000-maps_x[prev_wp];
    double center_y = 2000-maps_y[prev_wp];
    double centerToPos = distance(center_x,center_y,x_x,x_y);
    double centerToRef = distance(center_x,center_y,proj_x,proj_y);
    
    if(centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }
    
    // calculate s value
    double frenet_s = 0;
    for(int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
    }
    
    frenet_s += distance(0,0,proj_x,proj_y);
    
    return {frenet_s,frenet_d};
    
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int prev_wp = -1;
    
    while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
        prev_wp++;
    }
    
    int wp2 = (prev_wp+1)%maps_x.size();
    
    double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s-maps_s[prev_wp]);
    
    double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
    double seg_y = maps_y[prev_wp]+seg_s*sin(heading);
    
    double perp_heading = heading-pi()/2;
    
    double x = seg_x + d*cos(perp_heading);
    double y = seg_y + d*sin(perp_heading);
    
    return {x,y};
    
}

// cost functions
// All cost functions has a range of 0 to 1 so that it is easy to assign weights to each cost

// 1. Inefficiency cost, the slower the front car is, the higher the cost
// This function does not consider the car's current lane
// When the nearest frontal car in a particular lane drives at speed limit, or if there is no car at all
// The cost is 0.
double inefficiency_cost(double s_self, double target_speed, int lane_self, int intended_lane, const vector<vector<double>> &sensor_fusion ){
    
    double currentMinS = INFINITY;
    int currentIndex = -1;
    
    double total_cost = 0;
    
    // find the closest car in the intended lane
    for(int i=0; i<sensor_fusion.size(); i++){
        double s_other = sensor_fusion[i][5];
        double d_other = sensor_fusion[i][6];
        
        // if the car is in that lane
        if(d_other > (4*intended_lane) && d_other < (4+4*intended_lane)){
            double distance = s_other - s_self;
            // if the car is in front of us
            if(distance < currentMinS && distance > 0){
                currentMinS = distance;
                currentIndex = i;
            }
        }
    }
    
    // if there is no car detected, add 0 to total cost
    if(currentIndex == -1){
        
    }
    else{
    
        // get the speed of that car, calculate the cost
        double vx_other = sensor_fusion[currentIndex][3];
        double vy_other = sensor_fusion[currentIndex][4];
        
        double speed_other = sqrt(vx_other*vx_other + vy_other*vy_other);
        
        // if that car is driving over the target speed, just add 0 to total cost
        if(speed_other > target_speed){
            
        }
        else{
           total_cost += (target_speed - speed_other)/target_speed;
        }
    }
    
    
    return total_cost;
    
}

// 2. Distance cost
// Instead of just detecting the cars in front of us, also detect the cars that are behind us
// This makes it hard to just reuse the code for inefficiency cost
// get the car that is closest to our car
// get the distance of that car, if that car is more than the safety distance (30m), the cost is 0.
// In order to minimize lane changing and prevent side collision, the distance cost of the current lane is halved

double distance_cost(double s_self, int lane_self, int intended_lane, const vector<vector<double>> &sensor_fusion){
    
    double currentMinS = INFINITY;
    int currentIndex = -1;
    
    double total_cost = 0;
    
    // find the closest car in the intended lane
    for(int i=0; i<sensor_fusion.size(); i++){
        double s_other = sensor_fusion[i][5];
        double d_other = sensor_fusion[i][6];
        
        // if the car is in that lane
        if(d_other > (4*intended_lane) && d_other < (4+4*intended_lane)){
            double distance = abs(s_other - s_self);
            // if the car is the closest we have found, no matter it is in front or behind us
            if(distance < currentMinS){
                currentMinS = distance;
                currentIndex = i;
            }
        }
    }
    
    // if there is no car detected, the cost is 0
    if(currentIndex == -1){
        
    }

    // if the car is in our lane, the cost is halved
    if(intended_lane == lane_self){
        if(currentMinS <= 30)
            total_cost += (30 - currentMinS) / 60;

    }
    else{
        // only consider the car that is close enough
        if(currentMinS <= 30)
            total_cost += (30 - currentMinS) / 30;
    }
    

    return total_cost;

}

// 3. Future distance cost
// Cost based on the distance of the other car to our car 2 seconds later ( It takes 2 seconds to change lane)
// Both our car and the other car's location are projected into 2 seconds later, assuming constant s-axis velocity and no lane changing behavior
// Same as distance cost, the cost for center lane is halved
double future_distance_cost(double s_self, double v_self, int lane_self, int intended_lane, const vector<vector<double>> &sensor_fusion){
    
    double currentMinS = INFINITY;
    int currentIndex = -1;
    
    double total_cost = 0;
    
    // find the closest car in the intended lane
    for(int i=0; i<sensor_fusion.size(); i++){
        double s_other = sensor_fusion[i][5];
        double d_other = sensor_fusion[i][6];
        
        double vx_other = sensor_fusion[i][3];
        double vy_other = sensor_fusion[i][4];
        double speed_other = sqrt(vx_other*vx_other + vy_other*vy_other);
        
        // if the car is in that lane
        if(d_other > (4*intended_lane) && d_other < (4+4*intended_lane)){
            s_other = s_other + 2 * speed_other;
            s_self = s_self + 2 * v_self;
            double distance = abs(s_other - s_self);
            // if the car is the closest we have found, no matter it is in front or behind us
            if(distance < currentMinS){
                currentMinS = distance;
                currentIndex = i;
            }
        }
    }
    
    // if there is no car detected, the cost is 0
    if(currentIndex == -1){
        
    }
    
    // if the car is in our lane, the cost is halved
    if(intended_lane == lane_self){
        if(currentMinS <= 30)
            total_cost += (30 - currentMinS) / 60;
        
    }
    else{
        // only consider the car that is close enough
        if(currentMinS <= 30)
            total_cost += (30 - currentMinS) / 30;
    }
    

    return total_cost;
    
}

int main() {
    uWS::Hub h;
    
    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;
    
    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;
    
    ifstream in_map_(map_file_.c_str(), ifstream::in);
    
    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }
    
    // car state
    int lane = 1;
    double ref_velocity = 0;
    
    
    
    h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &ref_velocity](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                                                                                         uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length && length > 2 && data[0] == '4' && data[1] == '2') {
            
            auto s = hasData(data);
            
            if (s != "") {
                auto j = json::parse(s);
                
                string event = j[0].get<string>();
                
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    
                    // Main car's localization Data
                    double car_x = j[1]["x"];
                    double car_y = j[1]["y"];
                    double car_s = j[1]["s"];
                    double car_d = j[1]["d"];
                    double car_yaw = j[1]["yaw"];
                    double car_speed = j[1]["speed"];
                    
                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];
                    
                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    auto sensor_fusion = j[1]["sensor_fusion"];
                    
                    json msgJson;
                    
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;
                    
                    // the simulator provides us with the information for previous path
                    int prev_size = previous_path_x.size();
                    
                    
                    // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                    
                    double dist_inc = 0.5;
                    
                    // Stage 1: Sensor Fusion to detect other cars
                    
                    // Lanes in Frenet coordiate
                    // d = (|| 2 | 6 | 10 ||)
                    
                    // if there was previous data, set the car's position at the end_point_s
                    if(prev_size > 0){
                        car_s = end_path_s;
                    }
                    
                    bool too_close = false;
                    
                    // loop through all the cars in the sensor fusion
                    for(int i=0; i<sensor_fusion.size(); i++){
                    
                        // for each car in the sensor fusion, check if it is in our lane
                        // The data format for each car is: [ id, x, y, vx, vy, s, d].
                        double d_other = sensor_fusion[i][6];
                        double d_self = 4*lane+2;
                        
                        // if that is in our lane, check its future s position using the detected velocity * 0.02 (refresh rate). The next state should be after previous path, so this should be 0.02*prev_size*speed
                        if( (d_other > d_self - 2) && (d_other < d_self + 2)){
                            double vx = sensor_fusion[i][3];
                            double vy = sensor_fusion[i][4];
                            double speed = sqrt(vx*vx + vy*vy);
                            double s_current = sensor_fusion[i][5];
                            double s_future = s_current + speed*prev_size*0.02;
                            
                            // if the future position is with in our safety distance (30 meters front), flag for slow down
                            // This is independent to the cost functions. We must slow down first no matter
                            // what the cost function is
                            if((s_future > car_s) && (s_future - car_s < 30)){
                                too_close = true;
                                
                                
                                // lane change here
                                // use cost function to decide which lane to change, or not to change at all
                                // consider the cost function for all 3 lanes
                                int currentLane = -1;
                                double minCost = INFINITY;
                                vector<double> costs;
                                
                                for(int i=0; i<3; i++){
                                    
                                    // In most countries, lane change over 2 lanes is illegal
                                    // In addition, changing lanes over 2 lanes will cause an unacceptable
                                    // acceleration and jerk
                                    if( (lane == 0 && i != 2) || (lane == 2 && i != 0) || lane == 1 ){
                                   
                                        double cost_inefficiency = inefficiency_cost(car_s,49.5, lane, i, sensor_fusion);
                                        double cost_distance = distance_cost(car_s, lane, i, sensor_fusion);
                                        double cost_future_distance = future_distance_cost(car_s, ref_velocity, lane, i, sensor_fusion);
                                        
                                        cout << cost_inefficiency << endl;
                                        cout << cost_distance << endl;
                                        cout << cost_future_distance << endl;
                                        cout << "============" << endl;
                                        
                                        // The three costs have equal weight
                                        double total_cost = cost_inefficiency + cost_distance + cost_future_distance;
                                        if(total_cost < minCost){
                                            minCost = total_cost;
                                            currentLane = i;
                                            
                                        }
                                        costs.push_back(total_cost);
                                    }
                                
                                    
                                }
                                
                                cout << costs[0] << " " << costs[1] << " " << costs[2] << endl;
                                
                                lane = currentLane;
                                
                                
                            }
                        }
                    
                        
                    }
                    
                    
                    // Stage 2: Trajectory Generation using spline.h
                    
                    // control the speed
                    if(too_close){
                        ref_velocity -= 0.2;
                    }
                    else if(ref_velocity < 49.5){
                        ref_velocity += 0.2;
                    }
                    
                    // Create a list of waypoints in car's coordinate spacing in 30 meters
                    vector<double> ptsx;
                    vector<double> ptsy;
                    
                    // store the car's current state here so that we can use them to transform back
                    // from car's coordinate to world coordinate
                    double ref_x = car_x;
                    double ref_y = car_y;
                    // the simulator gives us this in degree, but we need radian
                    double ref_yaw = deg2rad(car_yaw);
                    
                    
                    
                    
                    
                    // if we do not have enough data from previous path, use car's current state
                    if(prev_size < 2)
                    {
                        // take the car's t-1 state
                        double x_prev = car_x - cos(ref_yaw);
                        double y_prev = car_y - sin(ref_yaw);
                        
                        // take the car's current state
                        
                        // push them to ptsx, ptsy
                        ptsx.push_back(x_prev);
                        ptsx.push_back(car_x);
                        ptsy.push_back(y_prev);
                        ptsy.push_back(car_y);
                        
                    }
                    // if we have enough data from previous path, use data from previous path
                    else
                    {
                        // take the car's state from previous path
                        // take the xy state at t-1
                        double x_t1 = previous_path_x[prev_size - 1];
                        double y_t1 = previous_path_y[prev_size - 1];
                        // take the xy state at t-2
                        double x_t2 = previous_path_x[prev_size - 2];
                        double y_t2 = previous_path_y[prev_size - 2];
                        // calculate the yaw between t-2 and t-1
                        double yaw_prev = atan2(y_t1 - y_t2 , x_t1 - x_t2);
                        // push t-1, t-2 into ptsx and ptsy
                        ptsx.push_back(x_t2);
                        ptsx.push_back(x_t1);
                        ptsy.push_back(y_t2);
                        ptsy.push_back(y_t1);
                        
                        // Because we used data from previous path, we need to set the ref_x and ref_y to the state at t-1, and ref_yaw between t-2 and t-1
                        ref_x = x_t1;
                        ref_y = y_t1;
                        ref_yaw = yaw_prev;
                        
                    }
                    
                    // take the car's 3 future states in Frenet in a 30m interval in s, convert it to XY, and push to the ptsx and ptsy
                    for(int i=30; i<=90; i+=30){
                        vector<double> f = getXY(car_s + i, 4*lane+2 , map_waypoints_s, map_waypoints_x, map_waypoints_y);
                        ptsx.push_back(f[0]);
                        ptsy.push_back(f[1]);
                    }
                
                    // convert all points in ptsx and ptsy from world coordinate to car's coordinate
                    // -Shift then -rotate
                    for(int i=0; i<ptsx.size(); i++){
                        double shift_x = ptsx[i] - ref_x;
                        double shift_y = ptsy[i] - ref_y;
                        
                        ptsx[i] = shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw);
                        ptsy[i] = shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw);
                        
                    }
                    
                    // fit the transformed ptsx and ptsy into a spline
                    tk::spline s;
                    s.set_points(ptsx, ptsy);
                    
                    
                    // push all the previous path to next state
                    for(int i=0; i<prev_size; i++){
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }
                    
                    // Determine the correct velocity
                    // The next waypoint is 30 meters away from us
                    double target_x = 30;
                    
                    // need to know the corresponding y value on the spline
                    double target_y = s(target_x);
                    // Euclidian distance
                    double target_dist = sqrt(target_x*target_x + target_y*target_y);
                    
                    
                    
                    
                    // the calculated X position so far
                    double x_add_on = 0;
                    
                    // Change the for loop to use previous results
                    for(int i = 0; i < 50 - prev_size; i++)
                    {
                        // Calculate N based on the reference velocity
                        double N = target_dist/(0.02 * ref_velocity/2.24);
                        
                        // Calculate X psoition based on previous position
                        // Get Y position using the spline
                        double x_point = x_add_on + (target_x)/N;
                        double y_point = s(x_point);
                        x_add_on = x_point;
                        
                        // Transform the points in car's coordinate back to world coordinate
                        // +Rotate then +shift
                        double x_ref = x_point;
                        double y_ref = y_point;
                        
                        x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
                        y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);
                        
                        x_point += ref_x;
                        y_point += ref_y;
                        
                        
                        
                        
                        // push XY coordinate in world coordinate into path planner
                        next_x_vals.push_back(x_point);
                        next_y_vals.push_back(y_point);
                    }
                    
                    
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;
                    
                    auto msg = "42[\"control\","+ msgJson.dump()+"]";
                    
                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                    
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });
    
    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });
    
    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });
    
    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
