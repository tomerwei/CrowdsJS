export function step(RADIUS, sceneEntities, world) {
    function distance(x1, y1, x2, y2) {
        return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }

    /*  -----------------------  */
    /*  TODO modify lines below  */
    /*  -----------------------  */
    function collisionConstraint(agent_i, agent_j) {
        const agentCentroidDist = distance(agent_i.px, agent_i.pz,
            agent_j.px, agent_j.pz);
        const agentDist = agentCentroidDist - AGENTSIZE;
        const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
        const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
        if (agentDist < 0) {
            sceneEntities[i].colliding = true;
            sceneEntities[j].colliding = true;
            agent_i.px += 0.5 * agentDist * dir_x
            agent_i.pz += 0.5 * agentDist * dir_z
            agent_j.px += -0.5 * agentDist * dir_x
            agent_j.pz += -0.5 * agentDist * dir_z
        }
    }
    /*  -----------------------  */

    const AGENTSIZE = RADIUS * 2;
    const epsilon = 0.0001;
    const timestep = 0.03;


    sceneEntities.forEach(function(item) {
        item.px = item.x + timestep * item.vx;
        item.pz = item.z + timestep * item.vz;
        item.py = item.y + timestep * item.vy;
    });


    var i = 0,
        j;
    while (i < sceneEntities.length) {
        j = i + 1;
        while (j < sceneEntities.length) {
            /*  TODO modify lines below if needed */
            /*  --------------------------------  */
            collisionConstraint(sceneEntities[i], sceneEntities[j])
            /*  --------------------------------  */
            j += 1;
        }
        i += 1
    }

    sceneEntities.forEach(function(item) {
        item.vx = (item.px - item.x) / timestep;
        item.vz = (item.pz - item.z) / timestep;
        item.vy = (item.py - item.y) / timestep;
        item.x = item.px;
        item.z = item.pz;
        item.y = item.py;

        if (item.x < -world.x / 2) {
            item.x = world.x / 2 - epsilon;
        } else if (item.x > world.x / 2) {
            item.x = -world.x / 2 + epsilon;
        }
        if (item.z < -world.z / 2) {
            item.z = world.z / 2 - epsilon;
        } else if (item.z > world.z / 2) {
            item.z = -world.z / 2 + epsilon;
        }
    });

}