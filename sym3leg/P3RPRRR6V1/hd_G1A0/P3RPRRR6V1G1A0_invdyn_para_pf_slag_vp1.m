% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:48
% EndTime: 2020-08-06 18:31:52
% DurationCPUTime: 4.17s
% Computational Cost: add. (16336->448), mult. (14877->662), div. (1575->13), fcn. (7791->86), ass. (0->305)
t216 = 2 * qJ(3,3);
t176 = sin(t216);
t157 = pkin(3) * t176;
t194 = sin(qJ(3,3));
t328 = 0.2e1 * t194;
t67 = 0.1e1 / (pkin(2) * t328 + t157 + (sin((pkin(7) + qJ(3,3))) + sin((-pkin(7) + qJ(3,3)))) * pkin(1));
t347 = t67 / 0.2e1;
t217 = 2 * qJ(3,2);
t177 = sin(t217);
t158 = pkin(3) * t177;
t196 = sin(qJ(3,2));
t327 = 0.2e1 * t196;
t68 = 0.1e1 / (pkin(2) * t327 + t158 + (sin((pkin(7) + qJ(3,2))) + sin((-pkin(7) + qJ(3,2)))) * pkin(1));
t346 = t68 / 0.2e1;
t218 = 2 * qJ(3,1);
t178 = sin(t218);
t159 = pkin(3) * t178;
t198 = sin(qJ(3,1));
t326 = 0.2e1 * t198;
t69 = 0.1e1 / (pkin(2) * t326 + t159 + (sin((pkin(7) + qJ(3,1))) + sin((-pkin(7) + qJ(3,1)))) * pkin(1));
t345 = t69 / 0.2e1;
t320 = Icges(3,2) / 0.2e1;
t344 = t320 - Icges(3,1) / 0.2e1;
t334 = -0.2e1 * pkin(2);
t261 = 0.2e1 * pkin(2);
t316 = m(3) / 0.2e1;
t175 = qJ(1,1) + pkin(7);
t190 = legFrame(1,3);
t138 = t190 + t175;
t120 = cos(t138);
t129 = t218 + t138;
t132 = -2 * qJ(3,1) + t138;
t156 = qJ(1,1) + t190;
t143 = qJ(3,1) + t156;
t144 = -qJ(3,1) + t156;
t262 = 0.2e1 * pkin(1);
t130 = qJ(3,1) + t138;
t131 = -qJ(3,1) + t138;
t265 = cos(t130) + cos(t131);
t268 = sin(t130) + sin(t131);
t213 = (-pkin(6) - pkin(5));
t325 = -2 * t213;
t54 = t265 * t261 + (cos(t144) + cos(t143)) * t262 + t268 * t325 + (cos(t132) + cos(t129) + 0.2e1 * t120) * pkin(3);
t342 = t54 * t345;
t174 = qJ(1,2) + pkin(7);
t189 = legFrame(2,3);
t137 = t189 + t174;
t119 = cos(t137);
t125 = t217 + t137;
t128 = -2 * qJ(3,2) + t137;
t155 = qJ(1,2) + t189;
t141 = qJ(3,2) + t155;
t142 = -qJ(3,2) + t155;
t126 = qJ(3,2) + t137;
t127 = -qJ(3,2) + t137;
t266 = cos(t126) + cos(t127);
t269 = sin(t126) + sin(t127);
t53 = t266 * t261 + (cos(t142) + cos(t141)) * t262 + t269 * t325 + (cos(t128) + cos(t125) + 0.2e1 * t119) * pkin(3);
t341 = t53 * t346;
t173 = qJ(1,3) + pkin(7);
t188 = legFrame(3,3);
t136 = t188 + t173;
t118 = cos(t136);
t121 = t216 + t136;
t124 = -2 * qJ(3,3) + t136;
t154 = qJ(1,3) + t188;
t139 = qJ(3,3) + t154;
t140 = -qJ(3,3) + t154;
t122 = qJ(3,3) + t136;
t123 = -qJ(3,3) + t136;
t267 = cos(t122) + cos(t123);
t270 = sin(t122) + sin(t123);
t52 = t267 * t261 + (cos(t140) + cos(t139)) * t262 + t270 * t325 + (cos(t124) + cos(t121) + 0.2e1 * t118) * pkin(3);
t340 = t52 * t347;
t117 = sin(t138);
t324 = 2 * t213;
t51 = t268 * t261 + (sin(t144) + sin(t143)) * t262 + t265 * t324 + (sin(t132) + sin(t129) + 0.2e1 * t117) * pkin(3);
t339 = t51 * t345;
t116 = sin(t137);
t50 = t269 * t261 + (sin(t142) + sin(t141)) * t262 + t266 * t324 + (sin(t128) + sin(t125) + 0.2e1 * t116) * pkin(3);
t338 = t50 * t346;
t115 = sin(t136);
t49 = t270 * t261 + (sin(t140) + sin(t139)) * t262 + t267 * t324 + (sin(t124) + sin(t121) + 0.2e1 * t115) * pkin(3);
t337 = t49 * t347;
t186 = sin(pkin(7));
t314 = pkin(1) * t186;
t336 = -t213 + t314;
t335 = -0.2e1 * pkin(1);
t211 = xDP(2);
t212 = xDP(1);
t222 = 0.1e1 / pkin(3);
t280 = t222 * t67;
t55 = t118 * t325 + sin(t154) * t335 + t115 * t334 - t270 * pkin(3);
t58 = t115 * t324 + cos(t154) * t335 + t118 * t334 - t267 * pkin(3);
t43 = (t211 * t55 + t212 * t58) * t280;
t40 = t43 ^ 2;
t279 = t222 * t68;
t56 = t119 * t325 + sin(t155) * t335 + t116 * t334 - t269 * pkin(3);
t59 = t116 * t324 + cos(t155) * t335 + t119 * t334 - t266 * pkin(3);
t44 = (t211 * t56 + t212 * t59) * t279;
t41 = t44 ^ 2;
t278 = t222 * t69;
t57 = t120 * t325 + sin(t156) * t335 + t117 * t334 - t268 * pkin(3);
t60 = t117 * t324 + cos(t156) * t335 + t120 * t334 - t265 * pkin(3);
t45 = (t211 * t57 + t212 * t60) * t278;
t42 = t45 ^ 2;
t187 = cos(pkin(7));
t160 = t187 * pkin(1);
t200 = cos(qJ(3,3));
t298 = t200 * pkin(3) + pkin(2);
t90 = t160 + t298;
t87 = 0.1e1 / t90;
t61 = (-t115 * t212 + t118 * t211) * t87;
t333 = t61 ^ 2;
t202 = cos(qJ(3,2));
t297 = t202 * pkin(3) + pkin(2);
t91 = t160 + t297;
t88 = 0.1e1 / t91;
t62 = (-t116 * t212 + t119 * t211) * t88;
t332 = t62 ^ 2;
t204 = cos(qJ(3,1));
t296 = t204 * pkin(3) + pkin(2);
t92 = t160 + t296;
t89 = 0.1e1 / t92;
t63 = (-t117 * t212 + t120 * t211) * t89;
t331 = t63 ^ 2;
t221 = pkin(3) ^ 2;
t224 = pkin(1) ^ 2;
t263 = pkin(2) ^ 2 + t224;
t239 = (t213 ^ 2) + t314 * t325 + t263;
t260 = 0.2e1 * t160;
t330 = t260 * t334 - t221 - 0.2e1 * t239;
t145 = t160 + pkin(2);
t329 = -0.2e1 * t145;
t323 = m(2) * g(3);
t322 = m(3) * rSges(3,2);
t321 = -Icges(3,5) / 0.4e1;
t319 = Icges(3,6) / 0.4e1;
t219 = rSges(3,2) ^ 2;
t220 = rSges(3,1) ^ 2;
t318 = (-t219 + t220) * t316 + t344;
t317 = t220 / 0.2e1;
t206 = rSges(3,3) + pkin(5);
t315 = m(3) * t222;
t195 = sin(qJ(1,3));
t313 = pkin(1) * t195;
t197 = sin(qJ(1,2));
t312 = pkin(1) * t197;
t199 = sin(qJ(1,1));
t311 = pkin(1) * t199;
t304 = t55 * t67;
t303 = t56 * t68;
t302 = t57 * t69;
t301 = t58 * t67;
t300 = t59 * t68;
t299 = t60 * t69;
t295 = t115 * t87;
t294 = t116 * t88;
t293 = t117 * t89;
t292 = t118 * t87;
t291 = t119 * t88;
t290 = t120 * t89;
t289 = t194 * t43;
t288 = t194 * t61;
t287 = t196 * t44;
t286 = t196 * t62;
t285 = t198 * t45;
t284 = t198 * t63;
t283 = t200 * t43;
t282 = t202 * t44;
t281 = t204 * t45;
t133 = (t219 + t220) * m(3) + Icges(3,3);
t277 = t133 * t222;
t276 = t194 * t145;
t275 = t196 * t145;
t274 = t198 * t145;
t273 = t206 * t186;
t215 = m(2) + m(3);
t272 = t215 / 0.2e1;
t271 = pkin(3) * t261;
t259 = rSges(3,1) * t322;
t264 = -t259 / 0.2e1 + Icges(3,4) / 0.2e1;
t95 = rSges(3,1) * t200 - rSges(3,2) * t194;
t258 = t95 * t316;
t97 = rSges(3,1) * t202 - rSges(3,2) * t196;
t257 = t97 * t316;
t99 = rSges(3,1) * t204 - rSges(3,2) * t198;
t256 = t99 * t316;
t255 = t95 * t315;
t254 = t97 * t315;
t253 = t99 * t315;
t252 = pkin(3) * t289;
t251 = pkin(3) * t287;
t250 = pkin(3) * t285;
t34 = (t55 * t255 + t49 * t272) * t67;
t35 = (t56 * t254 + t50 * t272) * t68;
t36 = (t57 * t253 + t51 * t272) * t69;
t249 = t36 + t35 + t34;
t37 = (t58 * t255 + t52 * t272) * t67;
t38 = (t59 * t254 + t53 * t272) * t68;
t39 = (t60 * t253 + t54 * t272) * t69;
t248 = t39 + t38 + t37;
t243 = t206 + t314;
t226 = -t243 * rSges(3,1) * m(3) + Icges(3,5);
t77 = -t243 * t322 + Icges(3,6);
t64 = t194 * t226 + t200 * t77;
t247 = t64 * t280;
t65 = t196 * t226 + t202 * t77;
t246 = t65 * t279;
t66 = t198 * t226 + t204 * t77;
t245 = t66 * t278;
t244 = t145 * t316;
t240 = rSges(3,1) * t244;
t238 = t160 / 0.2e1 + pkin(2) / 0.2e1;
t237 = t314 / 0.4e1 + t206 / 0.4e1;
t182 = t200 ^ 2;
t74 = 0.1e1 / (t157 + 0.2e1 * t276);
t16 = (-t336 * t252 + (t182 * t221 + t200 * t271 + t298 * t260 + t239) * t61) / (t157 / 0.2e1 + t276) * t61 * t222 + 0.2e1 * (t90 * t283 - t288 * t336) * t74 * t43;
t94 = rSges(3,1) * t194 + rSges(3,2) * t200;
t236 = -t16 * t95 - t40 * t94;
t183 = t202 ^ 2;
t75 = 0.1e1 / (t158 + 0.2e1 * t275);
t17 = (-t336 * t251 + (t183 * t221 + t202 * t271 + t297 * t260 + t239) * t62) / (t158 / 0.2e1 + t275) * t62 * t222 + 0.2e1 * (t91 * t282 - t286 * t336) * t75 * t44;
t96 = rSges(3,1) * t196 + rSges(3,2) * t202;
t235 = -t17 * t97 - t41 * t96;
t184 = t204 ^ 2;
t76 = 0.1e1 / (t159 + 0.2e1 * t274);
t18 = (-t336 * t250 + (t184 * t221 + t204 * t271 + t296 * t260 + t239) * t63) / (t159 / 0.2e1 + t274) * t63 * t222 + 0.2e1 * (t92 * t281 - t284 * t336) * t76 * t45;
t98 = rSges(3,1) * t198 + rSges(3,2) * t204;
t234 = -t18 * t99 - t42 * t98;
t233 = t61 * t238;
t232 = t62 * t238;
t231 = t63 * t238;
t230 = pkin(2) + t95;
t229 = pkin(2) + t97;
t228 = pkin(2) + t99;
t227 = t206 ^ 2 + t219 / 0.2e1 + t317 + t263;
t225 = Icges(1,3) + Icges(2,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(3,1) / 0.2e1 + t320 + (t224 + (rSges(2,1) * t187 - rSges(2,2) * t186) * t262 + rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t205 = cos(qJ(1,1));
t203 = cos(qJ(1,2));
t201 = cos(qJ(1,3));
t193 = xDDP(1);
t192 = xDDP(2);
t191 = xDDP(3);
t181 = cos(t218);
t180 = cos(t217);
t179 = cos(t216);
t172 = t205 * pkin(1);
t170 = t203 * pkin(1);
t168 = t201 * pkin(1);
t166 = cos(t190);
t165 = cos(t189);
t164 = cos(t188);
t163 = sin(t190);
t162 = sin(t189);
t161 = sin(t188);
t153 = cos(t175);
t152 = cos(t174);
t151 = cos(t173);
t150 = sin(t175);
t149 = sin(t174);
t148 = sin(t173);
t146 = -Icges(3,4) + t259;
t93 = rSges(3,2) * t244;
t86 = (t317 - t219 / 0.2e1) * m(3) + t344;
t85 = g(1) * t166 + g(2) * t163;
t84 = g(1) * t165 + g(2) * t162;
t83 = g(1) * t164 + g(2) * t161;
t82 = -g(1) * t163 + g(2) * t166;
t81 = -g(1) * t162 + g(2) * t165;
t80 = -g(1) * t161 + g(2) * t164;
t48 = t181 * t318 + (t99 * t261 + (t228 * t187 + t273) * t262 + t227) * m(3) - t146 * t178 + t225;
t47 = t180 * t318 + (t97 * t261 + (t229 * t187 + t273) * t262 + t227) * m(3) - t146 * t177 + t225;
t46 = t179 * t318 + (t95 * t261 + (t230 * t187 + t273) * t262 + t227) * m(3) - t146 * t176 + t225;
t33 = t57 * t245 + t48 * t290;
t32 = t56 * t246 + t47 * t291;
t31 = t55 * t247 + t46 * t292;
t30 = t60 * t245 - t48 * t293;
t29 = t59 * t246 - t47 * t294;
t28 = t58 * t247 - t46 * t295;
t27 = t66 * t290 + (t51 * t256 + t57 * t277) * t69;
t26 = t65 * t291 + (t50 * t257 + t56 * t277) * t68;
t25 = t64 * t292 + (t49 * t258 + t55 * t277) * t67;
t24 = -t66 * t293 + (t54 * t256 + t60 * t277) * t69;
t23 = -t65 * t294 + (t53 * t257 + t59 * t277) * t68;
t22 = -t64 * t295 + (t52 * t258 + t58 * t277) * t67;
t21 = (t336 * t63 - 0.2e1 * t250) * t89 * t63;
t20 = (t336 * t62 - 0.2e1 * t251) * t88 * t62;
t19 = (t336 * t61 - 0.2e1 * t252) * t87 * t61;
t15 = (t331 * t204 * t330 + (-0.2e1 * t42 * t92 + ((-t181 * t92 + t184 * t329 - t145) * t63 + (t204 * t326 + t178) * t45 * t336) * t63) * pkin(3)) * t76;
t14 = (t332 * t202 * t330 + (-0.2e1 * t41 * t91 + ((-t180 * t91 + t183 * t329 - t145) * t62 + (t202 * t327 + t177) * t44 * t336) * t62) * pkin(3)) * t75;
t13 = (t333 * t200 * t330 + (-0.2e1 * t40 * t90 + ((-t179 * t90 + t182 * t329 - t145) * t61 + (t200 * t328 + t176) * t43 * t336) * t61) * pkin(3)) * t74;
t12 = t215 * t15;
t11 = t215 * t14;
t10 = t215 * t13;
t9 = -t323 - t12 + (-g(3) + t234) * m(3);
t8 = -t323 - t11 + (-g(3) + t235) * m(3);
t7 = -t323 - t10 + (-g(3) + t236) * m(3);
t6 = -t66 * t21 - t133 * t18 + 0.2e1 * (t146 * t184 + (t198 * t86 + t93) * t204 + t198 * t240 + t264) * t331 + ((-g(3) - t15) * t99 + (t82 * t150 + t85 * t153) * t98) * m(3);
t5 = -t65 * t20 - t133 * t17 + 0.2e1 * (t146 * t183 + (t196 * t86 + t93) * t202 + t196 * t240 + t264) * t332 + ((-g(3) - t14) * t97 + (t81 * t149 + t84 * t152) * t96) * m(3);
t4 = -t64 * t19 - t133 * t16 + 0.2e1 * (t146 * t182 + (t194 * t86 + t93) * t200 + t194 * t240 + t264) * t333 + ((-t13 - g(3)) * t95 + (t80 * t148 + t83 * t151) * t94) * m(3);
t3 = -t48 * t21 - t66 * t18 - 0.4e1 * t45 * ((t284 * t318 + t45 * t321) * t204 + t285 * t319 + (t184 - 0.1e1 / 0.2e1) * t63 * t146 + ((t204 * t231 - t237 * t285) * rSges(3,2) + (t198 * t231 + t237 * t281) * rSges(3,1)) * m(3)) - m(1) * (t85 * (-rSges(1,1) * t199 - rSges(1,2) * t205) + t82 * (rSges(1,1) * t205 - rSges(1,2) * t199)) - m(2) * (t85 * (-rSges(2,1) * t150 - rSges(2,2) * t153 - t311) + t82 * (rSges(2,1) * t153 - rSges(2,2) * t150 + t172)) - m(3) * (-t85 * t311 + t82 * t172 + (t85 * t206 + t82 * t228) * t153 + (t82 * t206 - t85 * t228) * t150);
t2 = -t47 * t20 - t65 * t17 - 0.4e1 * t44 * ((t286 * t318 + t44 * t321) * t202 + t287 * t319 + (t183 - 0.1e1 / 0.2e1) * t62 * t146 + ((t202 * t232 - t237 * t287) * rSges(3,2) + (t196 * t232 + t237 * t282) * rSges(3,1)) * m(3)) - m(1) * (t84 * (-rSges(1,1) * t197 - rSges(1,2) * t203) + t81 * (rSges(1,1) * t203 - rSges(1,2) * t197)) - m(2) * (t84 * (-rSges(2,1) * t149 - rSges(2,2) * t152 - t312) + t81 * (rSges(2,1) * t152 - rSges(2,2) * t149 + t170)) - m(3) * (-t84 * t312 + t81 * t170 + (t84 * t206 + t81 * t229) * t152 + (t81 * t206 - t84 * t229) * t149);
t1 = -t46 * t19 - t64 * t16 - 0.4e1 * t43 * ((t288 * t318 + t43 * t321) * t200 + t289 * t319 + (t182 - 0.1e1 / 0.2e1) * t61 * t146 + ((t200 * t233 - t237 * t289) * rSges(3,2) + (t194 * t233 + t237 * t283) * rSges(3,1)) * m(3)) - m(1) * (t83 * (-rSges(1,1) * t195 - rSges(1,2) * t201) + t80 * (rSges(1,1) * t201 - rSges(1,2) * t195)) - m(2) * (t83 * (-rSges(2,1) * t148 - rSges(2,2) * t151 - t313) + t80 * (rSges(2,1) * t151 - rSges(2,2) * t148 + t168)) - m(3) * (-t83 * t313 + t80 * t168 + (t83 * t206 + t80 * t230) * t151 + (t80 * t206 - t83 * t230) * t148);
t70 = [-t1 * t295 - t2 * t294 - t3 * t293 - m(4) * g(1) + (t6 * t299 + t5 * t300 + t4 * t301) * t222 + t9 * t342 + t8 * t341 + t7 * t340 + (t28 * t292 + t29 * t291 + t30 * t290 + (t22 * t304 + t23 * t303 + t24 * t302) * t222 + t39 * t339 + t38 * t338 + t37 * t337) * t192 + t248 * t191 + (-t28 * t295 - t29 * t294 - t30 * t293 + m(4) + (t22 * t301 + t23 * t300 + t24 * t299) * t222 + t39 * t342 + t38 * t341 + t37 * t340) * t193; t1 * t292 + t2 * t291 + t3 * t290 - m(4) * g(2) + (t6 * t302 + t5 * t303 + t4 * t304) * t222 + t9 * t339 + t8 * t338 + t7 * t337 + (-t31 * t295 - t32 * t294 - t33 * t293 + (t25 * t301 + t26 * t300 + t27 * t299) * t222 + t36 * t342 + t35 * t341 + t34 * t340) * t193 + t249 * t191 + (t31 * t292 + t32 * t291 + t33 * t290 + m(4) + (t25 * t304 + t26 * t303 + t27 * t302) * t222 + t36 * t339 + t35 * t338 + t34 * t337) * t192; -t10 - t11 - t12 + (m(4) + 0.3e1 * t215) * t191 + (-0.3e1 * m(2) - m(4)) * g(3) + t248 * t193 + t249 * t192 + (-0.3e1 * g(3) + t234 + t235 + t236) * m(3);];
tauX  = t70;
