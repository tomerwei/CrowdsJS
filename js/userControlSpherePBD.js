import * as THREE from "three";

export function distance(x1, y1, x2, y2) {
    return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

export function step(RADIUS, sceneEntities, world, scene, customParams = {}) {
    const AGENTSIZE = RADIUS * 2;
    const epsilon = 0.0001;
    const timestep = 0.03;
    const ITERNUM = 1; // 3
    const agentLength = RADIUS;

    // collision functions

    function rotateLineSegment(x1, y1, x2, y2, r) {
        // Calculate the center of the line segment
        const centerX = (x1 + x2) / 2;
        const centerY = (y1 + y2) / 2;

        // Translate the line segment so that its center is at the origin
        const x1p = x1 - centerX;
        const y1p = y1 - centerY;
        const x2p = x2 - centerX;
        const y2p = y2 - centerY;

        // Rotate the line segment about the origin
        const cosR = Math.cos(r);
        const sinR = Math.sin(r);
        const x1r = x1p * cosR - y1p * sinR;
        const y1r = x1p * sinR + y1p * cosR;
        const x2r = x2p * cosR - y2p * sinR;
        const y2r = x2p * sinR + y2p * cosR;

        // Translate the line segment back to its original position
        const newX1 = x1r + centerX;
        const newY1 = y1r + centerY;
        const newX2 = x2r + centerX;
        const newY2 = y2r + centerY;

        // Return the new endpoints of the line segment
        return [newX1, newY1, newX2, newY2];
    }

    function ClosestPointOnLineSegment(A, B, Point) {
        const AB = B.clone().sub(A);
        const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
        return A.clone().add(
            AB.clone().multiplyScalar(Math.min(Math.max(t, 0), 1))
        );
    }

    function projectPointOnLineSegment(A, B, Point) {
        const AB = B.clone().sub(A);
        const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
        return A.clone().add(
            AB.clone().multiplyScalar(Math.min(Math.max(t, 0), 1))
        );
    }

    function is_colliding_torso(x11, y11, x12, y12, x21, y21, x22, y22) {
        // console.log(segments_distance(x11, y11, x12, y12, x21, y21, x22, y22));
        return (
            segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) < 2 * RADIUS
        );
    }

    function segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) {
        // distance between two segments in the plane:
        // one segment is (x11, y11) to (x12, y12)
        // the other is (x21, y21) to (x22, y22)

        if (segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22)) return 0;

        // try each of the 4 vertices w/the other segment
        let distances = [
            point_segment_distance(x11, y11, x21, y21, x22, y22),
            point_segment_distance(x12, y12, x21, y21, x22, y22),
            point_segment_distance(x21, y21, x11, y11, x12, y12),
            point_segment_distance(x22, y22, x11, y11, x12, y12),
        ];
        return Math.min(...distances);
    }

    function segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22) {
        // whether two segments in the plane intersect:
        // one segment is (x11, y11) to (x12, y12)
        // the other is (x21, y21) to (x22, y22)

        let dx1 = x12 - x11;
        let dy1 = y12 - y11;
        let dx2 = x22 - x21;
        let dy2 = y22 - y21;

        let delta = dx2 * dy1 - dy2 * dx1;
        if (delta === 0) return false; // parallel segments

        let s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta;
        let t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / -delta;

        return 0 <= s && s <= 1 && 0 <= t && t <= 1;
    }

    function point_segment_distance(px, py, x1, y1, x2, y2) {
        let dx = x2 - x1;
        let dy = y2 - y1;

        if (dx === 0 && dy === 0) {
            // the segment's just a point
            return Math.hypot(px - x1, py - y1);
        }

        // Calculate the t that minimizes the distance.
        let t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

        // See if this represents one of the segment's
        // end points or a point in the middle.
        if (t < 0) {
            dx = px - x1;
            dy = py - y1;
        } else if (t > 1) {
            dx = px - x2;
            dy = py - y2;
        } else {
            let near_x = x1 + t * dx;
            let near_y = y1 + t * dy;
            dx = px - near_x;
            dy = py - near_y;
        }
        return Math.hypot(dx, dy);
    }

    /*  -----------------------  */
    /*  TODO modify lines below  */
    /*  -----------------------  */
    function distanceConstraint(agent_i, agent_j, desiredDistance) {
        const agentCentroidDist = distance(
            agent_i.px,
            agent_i.pz,
            agent_j.px,
            agent_j.pz
        );
        const agentDist = agentCentroidDist - desiredDistance;
        const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
        const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
        const agent_i_scaler =
            ((0.1 * agent_i.invmass) / (agent_i.invmass + agent_j.invmass)) *
            agentDist;
        const agent_j_scaler =
            ((0.1 * agent_j.invmass) / (agent_i.invmass + agent_j.invmass)) *
            agentDist;
        if (Math.abs(agentDist) > epsilon) {
            agent_i.px += agent_i_scaler * dir_x;
            agent_i.pz += agent_i_scaler * dir_z;
            agent_j.px += -agent_j_scaler * dir_x;
            agent_j.pz += -agent_j_scaler * dir_z;
        }
    }

    function collisionConstraint(agent_i, agent_j) {
        const agentCentroidDist = distance(
            agent_i.px,
            agent_i.pz,
            agent_j.px,
            agent_j.pz
        );
        const agentDist = agentCentroidDist - AGENTSIZE;
        const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
        const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
        const agent_i_scaler =
            (agent_i.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
        const agent_j_scaler =
            (agent_j.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
        if (agentDist < 0) {
            agent_i.px += agent_i_scaler * dir_x;
            agent_i.pz += agent_i_scaler * dir_z;
            agent_j.px += -agent_j_scaler * dir_x;
            agent_j.pz += -agent_j_scaler * dir_z;
        }
    }

    function agentVelocityPlanner() {
        sceneEntities.forEach(function (agent_i) {
            const distToGoal = distance(
                agent_i.x,
                agent_i.z,
                agent_i.goal_x,
                agent_i.goal_z
            );
            if (distToGoal > RADIUS) {
                const dir_x = (agent_i.goal_x - agent_i.x) / distToGoal;
                const dir_z = (agent_i.goal_z - agent_i.z) / distToGoal;
                agent_i.vx = agent_i.v_pref * dir_x;
                agent_i.vz = agent_i.v_pref * dir_z;
            }
            agent_i.vx = 0.9999 * agent_i.vx;
            agent_i.vz = 0.9999 * agent_i.vz;
        });
    }



    function longRangeConstraint(agent_i, agent_j) {
        const agentCentroidDist = distance(agent_i.px, agent_i.pz,
            agent_j.px, agent_j.pz);
        const radius_init = 2 * AGENTSIZE;
        const radius_sq_init = radius_init * radius_init;
        var radius_sq = radius_sq_init;
        const dv_i = 1.;  // 1./delta_t;
        let delta_correction_i = {"x":0, "y":0};
        let delta_correction_j= {"x":0, "y":0};
        if (agentCentroidDist < radius_init) {
            radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
        }
        const v_x = (agent_i.px - agent_i.x) / timestep - (agent_j.px - agent_j.x) / timestep;
        const v_y = (agent_i.pz - agent_i.z) / timestep - (agent_j.pz - agent_j.z) / timestep;
        const x0 = agent_i.x - agent_j.x;
        const y0 = agent_i.z - agent_j.z;
        const v_sq = v_x * v_x + v_y * v_y;
        const x0_sq = x0 * x0;
        const y0_sq = y0 * y0;
        const x_sq = x0_sq + y0_sq;
        const a = v_sq;
        const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
        const b_sq = b * b;
        const c = x_sq - radius_sq;
        const d_sq = b_sq - a * c;
        const d = Math.sqrt(d_sq);
        const tao = (b - d) / a;




        // if (agent_i.index % 2 !== 0 && agent_j.index % 2 !== 0){
        //
        //   C_TAO0 = 6;
        // }else if (agent_i.index % 2 === 0 && agent_j.index % 2 === 0){
        //
        //   C_TAO0 = 80;
        // }else {
        //   return;
        // }

        let lengthV;
        if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {
            const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
            const tao_sq = c_tao * c_tao;
            const grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
            const grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
            const grad_x_j = -grad_x_i;
            const grad_y_j = -grad_y_i;
            const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
            const s = stiff * tao_sq / (agent_i.invmass * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + agent_j.invmass * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed

            lengthV = Math.sqrt(s * agent_i.invmass * grad_x_i * s * agent_i.invmass * grad_x_i
                + s * agent_i.invmass * grad_y_i * s * agent_i.invmass * grad_y_i);

            delta_correction_i = clamp2D(s * agent_i.invmass * grad_x_i,
                s * agent_i.invmass * grad_y_i,
                MAX_DELTA);

            delta_correction_j = clamp2D(s * agent_j.invmass * grad_x_j,
                s * agent_j.invmass * grad_y_j,
                MAX_DELTA);
            agent_i.px += delta_correction_i.x;
            agent_i.pz += delta_correction_i.y;
            agent_j.px += delta_correction_j.x;
            agent_j.pz += delta_correction_j.y;

            agent_i.grad[0] += delta_correction_i.x;
            agent_i.grad[1] += delta_correction_i.y;
            agent_j.grad[0] += delta_correction_j.x;
            agent_j.grad[1] += delta_correction_j.y;

        }
    }

    let C_TAU_MAX = 20;
    // const C_MAX_ACCELERATION = 0.01;
    let C_TAO0 = 25; //
    const C_LONG_RANGE_STIFF = 0.02;
    const MAX_DELTA = 0.9;

    function clamp2D(vx,vy, maxValue) {
        const lengthV = Math.sqrt(vx * vx + vy * vy);
        if (lengthV > maxValue) {
            const mult = (maxValue / lengthV);
            vx *= mult;
            vy *= mult;
        }
        return {"x":vx, "y":vy}
    }


    function longRangeConstraintCapsule(best_i, best_j, p_best_i, p_best_j) {
        const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

        const radius_init = 2 * AGENTSIZE;
        const radius_sq_init = radius_init * radius_init;
        let radius_sq = radius_sq_init;
        const dv_i = 1.;  // 1./delta_t;
        let delta_correction_i = {"x":0, "y":0};
        let delta_correction_j= {"x":0, "y":0};
        if (agentCentroidDist < radius_init) {
            radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
        }
        const v_x = (p_best_i.x - best_i.x) / timestep - (p_best_j.x - best_j.x) / timestep;
        const v_y = (p_best_i.z - best_i.z) / timestep - (p_best_j.z - best_j.z) / timestep;
        const x0 = best_i.x - best_j.x;
        const y0 = best_i.z - best_j.z;
        const v_sq = v_x * v_x + v_y * v_y;
        const x0_sq = x0 * x0;
        const y0_sq = y0 * y0;
        const x_sq = x0_sq + y0_sq;
        const a = v_sq;
        const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
        const b_sq = b * b;
        const c = x_sq - radius_sq;
        const d_sq = b_sq - a * c;
        const d = Math.sqrt(d_sq);
        const tao = (b - d) / a;
        // console.log("ttc in long range paper: " + tao);

        let lengthV;
        let grad_x_i;
        let grad_y_i;
        let grad_x_j;
        let grad_y_j;
        let s;

        if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {
            const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
            const tao_sq = c_tao * c_tao;

            grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
            grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
            grad_x_j = -grad_x_i;
            grad_y_j = -grad_y_i;

            const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
            s = stiff * tao_sq / (0.5 * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + 0.5 * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed
            // console.log()

            delta_correction_i = clamp2D(s * 0.5 * grad_x_i,
                s * 0.5 * grad_y_i,
                MAX_DELTA);

            delta_correction_j = clamp2D(s * 0.5 * grad_x_j,
                s * 0.5 * grad_y_j,
                MAX_DELTA);

            // console.log("Long Range active");

        }else {
            grad_x_i = 0;
            grad_y_i = 0;
            grad_x_j = 0;
            grad_y_j = 0;
            s=0;
        }


        // return tao;
        return [
            delta_correction_i,
            delta_correction_j,
            [grad_x_i, grad_y_i],
            [grad_x_j, grad_y_j],
            s
        ];
    }

    function getBestPointWithWall(xi, zi, wall){

        const iCoords = rotateLineSegment(
            xi,
            zi + agentLength + RADIUS,
            xi,
            zi - agentLength - RADIUS,
            sceneEntities[i].agent.rotation.z
        );



        // Agent A
        const a = {
            tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
            base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
            radius: RADIUS,
        };
        // Wall B
        const b = {
            tip: new THREE.Vector3(wall.tip.x, 0, wall.tip.z),
            base: new THREE.Vector3(wall.base.x, 0, wall.base.z),
            radius: RADIUS,
        };


        // capsule A:
        const a_Normal = a.tip.clone().sub(a.base.clone()).normalize();
        const a_LineEndOffset = a_Normal.clone().multiplyScalar(a.radius);
        const a_A = a.base.clone().add(a_LineEndOffset);
        const a_B = a.tip.clone().sub(a_LineEndOffset);

        // capsule B:
        const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
        const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
        const b_A = b.base.clone().add(b_LineEndOffset);
        const b_B = b.tip.clone().sub(b_LineEndOffset);

        // vectors between line endpoints:
        const v0 = b_A.clone().sub(a_A);
        const v1 = b_B.clone().sub(a_A);
        const v2 = b_A.clone().sub(a_B);
        const v3 = b_B.clone().sub(a_B);

        // squared distances:
        const d0 = v0.clone().dot(v0);
        const d1 = v1.clone().dot(v1);
        const d2 = v2.clone().dot(v2);
        const d3 = v3.clone().dot(v3);

        // select best potential endpoint on capsule A:
        let bestA;
        if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
            bestA = a_B;
        } else {
            bestA = a_A;
        }

        // select point on capsule B line segment nearest to best potential endpoint on A capsule:
        const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);

        // now do the same for capsule A segment:
        bestA = ClosestPointOnLineSegment(a_A, a_B, bestB);

        return [bestA, bestB, a, b]
    }

    function getCircleCenterWithWall(xi, zi, wall){



        // Agent A
        const a = {
            tip: new THREE.Vector3(xi, 0, zi),
            base: new THREE.Vector3(xi, 0, zi),
            radius: RADIUS,
        };
        // Wall B
        const b = {
            tip: new THREE.Vector3(wall.tip.x, 0, wall.tip.z),
            base: new THREE.Vector3(wall.base.x, 0, wall.base.z),
            radius: RADIUS,
        };


        // capsule A:
        const a_A = a.base;
        const a_B = a.tip;

        // capsule B:
        const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
        const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
        const b_A = b.base.clone().add(b_LineEndOffset);
        const b_B = b.tip.clone().sub(b_LineEndOffset);

        // vectors between line endpoints:
        const v0 = b_A.clone().sub(a_A);
        const v1 = b_B.clone().sub(a_A);
        const v2 = b_A.clone().sub(a_B);
        const v3 = b_B.clone().sub(a_B);

        // squared distances:
        const d0 = v0.clone().dot(v0);
        const d1 = v1.clone().dot(v1);
        const d2 = v2.clone().dot(v2);
        const d3 = v3.clone().dot(v3);

        // select best potential endpoint on capsule A:
        let bestA;
        if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
            bestA = a_B;
        } else {
            bestA = a_A;
        }

        // this line will be replaced by projecting point to line segment:
        const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);



        return [bestA, bestB, a, b]
    }


    function dotProduct(vector1, vector2) {
        let result = 0;
        for (let i = 0; i < vector1.length; i++) {
            result += vector1[i] * vector2[i];
        }
        return result;
    }

    function areCollinear(vector1, vector2) {
        // Ensure vectors are of the same dimension
        if (vector1.length !== vector2.length) {
            return false;
        }

        // Find the ratio of the first non-zero pair of elements
        let ratio;
        for (let i = 0; i < vector1.length; i++) {
            if (vector1[i] !== 0 && vector2[i] !== 0) {
                ratio = vector1[i] / vector2[i];
                break;
            }
        }

        // Check if all corresponding elements are in the same ratio
        for (let i = 0; i < vector1.length; i++) {
            // Handle division by zero cases
            if (vector1[i] === 0 && vector2[i] !== 0 || vector1[i] !== 0 && vector2[i] === 0) {
                return false;
            }

            // Check the ratio
            if (vector1[i] !== 0 && vector2[i] !== 0) {
                if (vector1[i] / vector2[i] !== ratio) {
                    return false;
                }
            }
        }

        return true;
    }

    function signNoP(n){
        if (n >= 0){
            return 1;
        }else {
            return -1;
        }

    }




    /*  -----------------------  */

    agentVelocityPlanner();

    sceneEntities.forEach(function (item) {
        item.px = item.x + timestep * item.vx;
        item.pz = item.z + timestep * item.vz;
        item.py = item.y + timestep * item.vy;
    });






    let pbdIters = 0;
    let isColliding;
    var agent_a,
        agent_b,
        desDistance,
        i,
        j,
        idx = 0;

    while (pbdIters < ITERNUM) {

        // wall collision (based on short range)
        i=0;
        while(i<sceneEntities.length)
        {
            j=0;
            while(j<customParams.wallData.length)
            {
                // let [p_bestA, w_bestB, p_agent_i,p_agent_j] = getBestPointWithWall(sceneEntities[i].px, sceneEntities[i].pz, customParams.wallData[j]);
                let [p_bestA, w_bestB, p_agent_i,p_agent_j] = getCircleCenterWithWall(sceneEntities[i].px, sceneEntities[i].pz, customParams.wallData[j]);

                let penetration_normal = p_bestA.clone().sub(w_bestB);
                const len = penetration_normal.length();
                penetration_normal.divideScalar(len); // normalize
                const penetration_depth = sceneEntities[i].radius + sceneEntities[i].radius - len;
                const intersects = penetration_depth > 0;
                if (intersects) {
                    sceneEntities[i].colliding = true;

                    sceneEntities[i].px += penetration_normal.x * 1 * penetration_depth;
                    sceneEntities[i].pz += penetration_normal.z * 1 * penetration_depth;

                }

                j+=1;
            }
            i+=1
        }

        i = 0;
        while (i < sceneEntities.length) {
            j = i + 1;
            while (j < sceneEntities.length) {

                // spheres
                let bestA, bestB;
                longRangeConstraint(sceneEntities[i], sceneEntities[j]);
                bestA = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
                bestB = new THREE.Vector3(sceneEntities[j].x, 0, sceneEntities[j].z);

                // utilities
                customParams.best[i][j] = [bestA, bestB]
                customParams.best[j][i] = [bestB, bestA]

                // short range collision
                // didn't correct position in real time
                // let penetration_normal = bestA.clone().sub(bestB);
                // const len = penetration_normal.length();
                // penetration_normal.divideScalar(len); // normalize
                // const penetration_depth = sceneEntities[i].radius + sceneEntities[j].radius - len;
                // const intersects = penetration_depth > 0;
                // if (intersects) {
                //   sceneEntities[i].colliding = true;
                //   sceneEntities[j].colliding = true;
                //
                //   sceneEntities[i].px += penetration_normal.x * 0.5 * penetration_depth;
                //   sceneEntities[i].pz += penetration_normal.y * 0.5 * penetration_depth;
                //
                //   sceneEntities[j].px +=
                //       -1 * penetration_normal.x * 0.5 * penetration_depth;
                //   sceneEntities[j].pz +=
                //       -1 * penetration_normal.y * 0.5 * penetration_depth;
                // }


                j += 1;
            }
            i += 1;
        }





        pbdIters += 1;
    }


    sceneEntities.forEach(function (item) {
        item.vx = (item.px - item.x) / timestep;
        item.vz = (item.pz - item.z) / timestep;
        item.vy = (item.py - item.y) / timestep;
        item.x = item.px;
        item.z = item.pz;
        item.y = item.py;
    });
}
