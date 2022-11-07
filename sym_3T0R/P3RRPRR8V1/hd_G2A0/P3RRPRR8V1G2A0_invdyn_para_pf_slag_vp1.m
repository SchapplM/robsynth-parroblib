% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:09
% EndTime: 2022-11-04 17:04:14
% DurationCPUTime: 4.57s
% Computational Cost: add. (17061->456), mult. (28683->734), div. (3201->12), fcn. (21807->44), ass. (0->304)
t193 = pkin(4) + qJ(3,3);
t179 = 0.1e1 / t193;
t220 = xDP(3);
t221 = xDP(2);
t222 = xDP(1);
t196 = legFrame(3,2);
t158 = sin(t196);
t161 = cos(t196);
t203 = sin(qJ(1,3));
t209 = cos(qJ(1,3));
t189 = cos(pkin(5));
t137 = pkin(2) * t189 + pkin(1);
t208 = cos(qJ(2,3));
t188 = sin(pkin(5));
t202 = sin(qJ(2,3));
t279 = t188 * t202;
t236 = pkin(2) * t279 - t137 * t208;
t83 = -t193 * t209 - t203 * t236;
t319 = pkin(2) * t188;
t95 = t137 * t202 + t208 * t319;
t55 = t158 * t95 + t161 * t83;
t56 = -t158 * t83 + t161 * t95;
t174 = qJ(2,3) + pkin(5);
t154 = cos(t174);
t168 = t208 * pkin(1);
t346 = pkin(2) * t154 + t168;
t86 = t193 * t203 + t209 * t346;
t40 = (t220 * t86 + t221 * t56 + t222 * t55) * t179;
t359 = 0.2e1 * t40;
t194 = pkin(4) + qJ(3,2);
t180 = 0.1e1 / t194;
t197 = legFrame(2,2);
t159 = sin(t197);
t162 = cos(t197);
t205 = sin(qJ(1,2));
t211 = cos(qJ(1,2));
t210 = cos(qJ(2,2));
t204 = sin(qJ(2,2));
t278 = t188 * t204;
t235 = pkin(2) * t278 - t137 * t210;
t84 = -t194 * t211 - t205 * t235;
t96 = t137 * t204 + t210 * t319;
t57 = t159 * t96 + t162 * t84;
t58 = -t159 * t84 + t162 * t96;
t176 = qJ(2,2) + pkin(5);
t156 = cos(t176);
t169 = t210 * pkin(1);
t345 = pkin(2) * t156 + t169;
t87 = t194 * t205 + t211 * t345;
t41 = (t220 * t87 + t221 * t58 + t222 * t57) * t180;
t358 = 0.2e1 * t41;
t195 = pkin(4) + qJ(3,1);
t181 = 0.1e1 / t195;
t198 = legFrame(1,2);
t160 = sin(t198);
t163 = cos(t198);
t207 = sin(qJ(1,1));
t213 = cos(qJ(1,1));
t212 = cos(qJ(2,1));
t206 = sin(qJ(2,1));
t277 = t188 * t206;
t234 = pkin(2) * t277 - t137 * t212;
t85 = -t195 * t213 - t207 * t234;
t97 = t137 * t206 + t212 * t319;
t59 = t160 * t97 + t163 * t85;
t60 = -t160 * t85 + t163 * t97;
t177 = qJ(2,1) + pkin(5);
t157 = cos(t177);
t170 = t212 * pkin(1);
t344 = pkin(2) * t157 + t170;
t88 = t195 * t207 + t213 * t344;
t42 = (t220 * t88 + t221 * t60 + t222 * t59) * t181;
t357 = 0.2e1 * t42;
t270 = 0.2e1 * pkin(1);
t269 = t189 * t270;
t231 = pkin(2) ^ 2;
t232 = pkin(1) ^ 2;
t271 = -t231 - t232;
t116 = pkin(2) * t269 - t271;
t106 = 0.1e1 / t346;
t73 = (t158 * t222 + t161 * t221) * t106;
t70 = t73 ^ 2;
t355 = t116 * t70;
t107 = 0.1e1 / t345;
t74 = (t159 * t222 + t162 * t221) * t107;
t71 = t74 ^ 2;
t354 = t116 * t71;
t108 = 0.1e1 / t344;
t75 = (t160 * t222 + t163 * t221) * t108;
t72 = t75 ^ 2;
t353 = t116 * t72;
t148 = sin(t174);
t165 = t202 * pkin(1);
t352 = rSges(3,1) * t148 + t154 * rSges(3,2) + t165;
t351 = pkin(2) * t148 + t165;
t150 = sin(t176);
t166 = t204 * pkin(1);
t350 = rSges(3,1) * t150 + rSges(3,2) * t156 + t166;
t349 = pkin(2) * t150 + t166;
t151 = sin(t177);
t167 = t206 * pkin(1);
t348 = rSges(3,1) * t151 + rSges(3,2) * t157 + t167;
t347 = pkin(2) * t151 + t167;
t337 = m(3) * rSges(3,2);
t267 = t188 * t337;
t98 = -m(2) * rSges(2,1) + t267 + (-rSges(3,1) * t189 - pkin(1)) * m(3);
t339 = m(2) * rSges(2,2);
t99 = t339 + (rSges(3,1) * t188 + rSges(3,2) * t189) * m(3);
t343 = -t99 * t206 - t212 * t98;
t342 = -t99 * t204 - t210 * t98;
t341 = -t99 * t202 - t208 * t98;
t340 = m(1) * rSges(1,1);
t338 = m(2) * rSges(2,3);
t89 = rSges(3,1) * t154 - rSges(3,2) * t148 + t168;
t336 = m(3) * t89;
t90 = rSges(3,1) * t156 - rSges(3,2) * t150 + t169;
t335 = m(3) * t90;
t91 = rSges(3,1) * t157 - rSges(3,2) * t151 + t170;
t334 = m(3) * t91;
t228 = rSges(2,2) ^ 2;
t230 = rSges(2,1) ^ 2;
t109 = t232 * m(3) + (-t228 + t230) * m(2) + Icges(2,2) - Icges(2,1);
t333 = t109 / 0.2e1;
t227 = rSges(3,2) ^ 2;
t229 = rSges(3,1) ^ 2;
t118 = m(3) * (-t227 + t229) - Icges(3,1) + Icges(3,2);
t332 = t118 / 0.2e1;
t331 = m(3) * t179;
t330 = m(3) * t180;
t329 = m(3) * t181;
t190 = rSges(3,3) + qJ(3,3);
t328 = m(3) * t190;
t191 = rSges(3,3) + qJ(3,2);
t327 = m(3) * t191;
t192 = rSges(3,3) + qJ(3,1);
t326 = m(3) * t192;
t266 = t203 * t319;
t288 = t137 * t203;
t61 = (-t158 * t288 + t161 * t319) * t208 + (t137 * t161 + t158 * t266) * t202;
t92 = 0.1e1 / t236;
t318 = t61 * t92;
t265 = t205 * t319;
t287 = t137 * t205;
t62 = (-t159 * t287 + t162 * t319) * t210 + (t137 * t162 + t159 * t265) * t204;
t93 = 0.1e1 / t235;
t317 = t62 * t93;
t264 = t207 * t319;
t286 = t137 * t207;
t63 = (-t160 * t286 + t163 * t319) * t212 + (t137 * t163 + t160 * t264) * t206;
t94 = 0.1e1 / t234;
t316 = t63 * t94;
t64 = (t158 * t319 + t161 * t288) * t208 + (t137 * t158 - t161 * t266) * t202;
t315 = t64 * t92;
t65 = (t159 * t319 + t162 * t287) * t210 + (t137 * t159 - t162 * t265) * t204;
t314 = t65 * t93;
t66 = (t160 * t319 + t163 * t286) * t212 + (t137 * t160 - t163 * t264) * t206;
t313 = t66 * t94;
t312 = t89 * t92;
t311 = t90 * t93;
t310 = t91 * t94;
t31 = (t209 * t220 - (t221 * t61 + t222 * t64) * t92) * t179;
t306 = t31 / 0.2e1;
t32 = (t211 * t220 - (t221 * t62 + t222 * t65) * t93) * t180;
t305 = t32 / 0.2e1;
t33 = (t213 * t220 - (t221 * t63 + t222 * t66) * t94) * t181;
t304 = t33 / 0.2e1;
t300 = t106 * t351;
t299 = t106 * t158;
t298 = t106 * t161;
t297 = t107 * t349;
t296 = t107 * t159;
t295 = t107 * t162;
t294 = t108 * t347;
t293 = t108 * t160;
t292 = t108 * t163;
t224 = 0.2e1 * qJ(2,3);
t182 = sin(t224);
t291 = t109 * t182;
t225 = 0.2e1 * qJ(2,2);
t183 = sin(t225);
t290 = t109 * t183;
t226 = 0.2e1 * qJ(2,1);
t184 = sin(t226);
t289 = t109 * t184;
t138 = 0.2e1 * t174;
t124 = cos(t138);
t142 = -rSges(3,1) * t337 + Icges(3,4);
t285 = t142 * t124;
t139 = 0.2e1 * t176;
t125 = cos(t139);
t284 = t142 * t125;
t140 = 0.2e1 * t177;
t126 = cos(t140);
t283 = t142 * t126;
t144 = rSges(2,1) * t339 - Icges(2,4);
t185 = cos(t224);
t282 = t144 * t185;
t186 = cos(t225);
t281 = t144 * t186;
t187 = cos(t226);
t280 = t144 * t187;
t276 = -rSges(2,1) * t338 + Icges(2,5);
t172 = pkin(5) + t224;
t152 = cos(t172);
t275 = t152 + t189;
t173 = pkin(5) + t226;
t153 = cos(t173);
t274 = t153 + t189;
t175 = t225 + pkin(5);
t155 = cos(t175);
t273 = t155 + t189;
t272 = t228 + t230;
t127 = rSges(3,2) * t328 - Icges(3,6);
t130 = rSges(3,1) * t328 - Icges(3,5);
t143 = rSges(2,2) * t338 - Icges(2,6);
t239 = pkin(1) * t328 - t276;
t46 = (t127 * t188 - t130 * t189 - t239) * t202 - (t127 * t189 + t130 * t188 + t143) * t208;
t263 = t179 * t46 * t92;
t128 = rSges(3,2) * t327 - Icges(3,6);
t131 = rSges(3,1) * t327 - Icges(3,5);
t238 = pkin(1) * t327 - t276;
t47 = (t128 * t188 - t131 * t189 - t238) * t204 - (t128 * t189 + t131 * t188 + t143) * t210;
t262 = t180 * t47 * t93;
t129 = rSges(3,2) * t326 - Icges(3,6);
t132 = rSges(3,1) * t326 - Icges(3,5);
t237 = pkin(1) * t326 - t276;
t48 = (t129 * t188 - t132 * t189 - t237) * t206 - (t129 * t189 + t132 * t188 + t143) * t212;
t261 = t181 * t48 * t94;
t260 = t106 * t209 * t46;
t259 = t107 * t211 * t47;
t258 = t108 * t213 * t48;
t257 = m(3) * t270;
t256 = pkin(2) * t270;
t249 = t227 / 0.2e1 + t229 / 0.2e1 + t232 / 0.2e1;
t146 = sin(t172);
t248 = rSges(3,1) * t146 + rSges(3,2) * t152;
t147 = sin(t173);
t247 = rSges(3,1) * t147 + rSges(3,2) * t153;
t149 = sin(t175);
t245 = rSges(3,1) * t149 + rSges(3,2) * t155;
t103 = g(1) * t161 - g(2) * t158;
t242 = g(3) * t209 + t103 * t203;
t104 = g(1) * t162 - g(2) * t159;
t241 = g(3) * t211 + t104 * t205;
t105 = g(1) * t163 - g(2) * t160;
t240 = g(3) * t213 + t105 * t207;
t136 = pkin(1) * t267;
t233 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 - t136 + (0.2e1 * rSges(2,3) ^ 2 + t272) * m(2) / 0.2e1;
t201 = xDDP(1);
t200 = xDDP(2);
t199 = xDDP(3);
t178 = g(3) * t340;
t141 = m(1) * rSges(1,2) - t338;
t133 = t141 * g(3);
t123 = sin(t140);
t122 = sin(t139);
t121 = sin(t138);
t102 = g(1) * t160 + g(2) * t163;
t101 = g(1) * t159 + g(2) * t162;
t100 = g(1) * t158 + g(2) * t161;
t79 = -0.2e1 * t136 + t272 * m(2) + Icges(2,3) + Icges(3,3) + (rSges(3,1) * t269 + t227 + t229 + t232) * m(3);
t69 = (-t213 * t91 + t88) * t329;
t68 = (-t211 * t90 + t87) * t330;
t67 = (-t209 * t89 + t86) * t331;
t45 = t126 * t332 + t187 * t333 + (t192 ^ 2 + (rSges(3,1) * t274 - rSges(3,2) * t147) * pkin(1) + t249) * m(3) + t142 * t123 - t144 * t184 + t233;
t44 = t125 * t332 + t186 * t333 + (t191 ^ 2 + (rSges(3,1) * t273 - rSges(3,2) * t149) * pkin(1) + t249) * m(3) + t142 * t122 - t144 * t183 + t233;
t43 = t124 * t332 + t185 * t333 + (t190 ^ 2 + (rSges(3,1) * t275 - rSges(3,2) * t146) * pkin(1) + t249) * m(3) + t142 * t121 - t144 * t182 + t233;
t39 = (t310 * t66 + t59) * t329;
t38 = (t311 * t65 + t57) * t330;
t37 = (t312 * t64 + t55) * t331;
t36 = (t310 * t63 + t60) * t329;
t35 = (t311 * t62 + t58) * t330;
t34 = (t312 * t61 + t56) * t331;
t30 = (t213 * t45 - t334 * t88) * t181;
t29 = (t211 * t44 - t335 * t87) * t180;
t28 = (t209 * t43 - t336 * t86) * t179;
t27 = -t261 * t66 + t293 * t79;
t26 = -t262 * t65 + t296 * t79;
t25 = -t263 * t64 + t299 * t79;
t24 = -t261 * t63 + t292 * t79;
t23 = -t262 * t62 + t295 * t79;
t22 = -t263 * t61 + t298 * t79;
t21 = t48 * t293 + (-t313 * t45 - t334 * t59) * t181;
t20 = t47 * t296 + (-t314 * t44 - t335 * t57) * t180;
t19 = t46 * t299 + (-t315 * t43 - t336 * t55) * t179;
t18 = t48 * t292 + (-t316 * t45 - t334 * t60) * t181;
t17 = t47 * t295 + (-t317 * t44 - t335 * t58) * t180;
t16 = t46 * t298 + (-t318 * t43 - t336 * t56) * t179;
t15 = (-0.1e1 / (t170 + (t189 * t212 - t277) * pkin(2)) * t353 + (t234 * t33 + t357) * t33) * t181;
t14 = (-0.1e1 / (t169 + (t189 * t210 - t278) * pkin(2)) * t354 + (t235 * t32 + t358) * t32) * t180;
t13 = (0.1e1 / (-t168 + (-t189 * t208 + t279) * pkin(2)) * t355 + (t236 * t31 + t359) * t31) * t179;
t12 = (-t353 + ((-t231 * t126 - t232 * t187 - t256 * t274 + t271) * t304 + t344 * t357 + 0.2e1 * (-t195 * t304 + t347 * t75) * t195) * t33) * t181;
t11 = (-t354 + ((-t231 * t125 - t232 * t186 - t256 * t273 + t271) * t305 + t345 * t358 + 0.2e1 * (-t194 * t305 + t349 * t74) * t194) * t32) * t180;
t10 = (-t355 + ((-t231 * t124 - t232 * t185 - t256 * t275 + t271) * t306 + t346 * t359 + 0.2e1 * (-t193 * t306 + t351 * t73) * t193) * t31) * t179;
t9 = -t48 * t15 + t79 * t72 * t294 + (t102 * t99 - t240 * t98) * t206 + t212 * (t102 * t98 + t240 * t99) + ((t123 * t332 - t283 + t289 / 0.2e1 + t280) * t33 + (t247 * t33 * pkin(1) - 0.2e1 * t42 * t348) * m(3)) * t33;
t8 = -t47 * t14 + t79 * t71 * t297 + (t101 * t99 - t241 * t98) * t204 + t210 * (t101 * t98 + t241 * t99) + ((t122 * t332 - t284 + t290 / 0.2e1 + t281) * t32 + (t245 * t32 * pkin(1) - 0.2e1 * t41 * t350) * m(3)) * t32;
t7 = -t46 * t13 + t79 * t70 * t300 + (t100 * t99 - t242 * t98) * t202 + t208 * (t100 * t98 + t242 * t99) + ((t121 * t332 - t285 + t291 / 0.2e1 + t282) * t31 + (t248 * t31 * pkin(1) - 0.2e1 * t40 * t352) * m(3)) * t31;
t6 = (-g(3) * t207 + t105 * t213 + t91 * t15 - t12 + (-t33 * t192 + 0.2e1 * t348 * t75) * t33) * m(3);
t5 = (-g(3) * t205 + t104 * t211 + t90 * t14 - t11 + (-t32 * t191 + 0.2e1 * t350 * t74) * t32) * m(3);
t4 = (-g(3) * t203 + t103 * t209 + t89 * t13 - t10 + (-t31 * t190 + 0.2e1 * t352 * t73) * t31) * m(3);
t3 = -t45 * t15 + t12 * t334 + (-g(3) * t326 + t133) * t213 + t207 * (t343 * g(3) + t178) + ((-t340 - t343) * t213 + t207 * (t141 - t326)) * t105 + (t129 * t151 - t132 * t157 + t143 * t206 - t212 * t237 + t294 * t48) * t72 + (t326 * t357 + (-t118 * t123 - t247 * t257 - 0.2e1 * t280 + 0.2e1 * t283 - t289) * t75) * t33;
t2 = -t44 * t14 + t11 * t335 + (-g(3) * t327 + t133) * t211 + t205 * (t342 * g(3) + t178) + ((-t340 - t342) * t211 + t205 * (t141 - t327)) * t104 + (t128 * t150 - t131 * t156 + t143 * t204 - t210 * t238 + t297 * t47) * t71 + (t327 * t358 + (-t118 * t122 - t245 * t257 - 0.2e1 * t281 + 0.2e1 * t284 - t290) * t74) * t32;
t1 = -t43 * t13 + t10 * t336 + (-g(3) * t328 + t133) * t209 + t203 * (t341 * g(3) + t178) + ((-t340 - t341) * t209 + t203 * (t141 - t328)) * t103 + (t127 * t148 - t130 * t154 + t143 * t202 - t208 * t239 + t300 * t46) * t70 + (t328 * t359 + (-t118 * t121 - t248 * t257 - 0.2e1 * t282 + 0.2e1 * t285 - t291) * t73) * t31;
t49 = [t7 * t299 + t8 * t296 + t9 * t293 - m(4) * g(1) + (t25 * t298 + t26 * t295 + t27 * t292) * t200 + (t25 * t299 + t26 * t296 + t27 * t293 + m(4)) * t201 + ((-t21 * t313 + t39 * t59) * t201 + (-t21 * t316 + t39 * t60) * t200 + (t21 * t213 + t39 * t88) * t199 - t3 * t313 + t59 * t6) * t181 + ((-t20 * t314 + t38 * t57) * t201 + (-t20 * t317 + t38 * t58) * t200 + (t20 * t211 + t38 * t87) * t199 - t2 * t314 + t57 * t5) * t180 + ((-t19 * t315 + t37 * t55) * t201 + (-t19 * t318 + t37 * t56) * t200 + (t19 * t209 + t37 * t86) * t199 - t1 * t315 + t55 * t4) * t179; t7 * t298 + t8 * t295 + t9 * t292 - m(4) * g(2) + (t22 * t299 + t23 * t296 + t24 * t293) * t201 + (t22 * t298 + t23 * t295 + t24 * t292 + m(4)) * t200 + ((-t18 * t313 + t36 * t59) * t201 + (-t18 * t316 + t36 * t60) * t200 + (t18 * t213 + t36 * t88) * t199 - t3 * t316 + t60 * t6) * t181 + ((-t17 * t314 + t35 * t57) * t201 + (-t17 * t317 + t35 * t58) * t200 + (t17 * t211 + t35 * t87) * t199 - t2 * t317 + t58 * t5) * t180 + ((-t16 * t315 + t34 * t55) * t201 + (-t16 * t318 + t34 * t56) * t200 + (t16 * t209 + t34 * t86) * t199 - t1 * t318 + t56 * t4) * t179; (-g(3) + t199) * m(4) + ((t160 * t258 - t30 * t313 + t59 * t69) * t201 + (t163 * t258 - t30 * t316 + t60 * t69) * t200 + (t213 * t30 + t69 * t88) * t199 + t213 * t3 + t88 * t6) * t181 + ((t159 * t259 - t29 * t314 + t57 * t68) * t201 + (t162 * t259 - t29 * t317 + t58 * t68) * t200 + (t211 * t29 + t68 * t87) * t199 + t211 * t2 + t87 * t5) * t180 + ((t158 * t260 - t28 * t315 + t55 * t67) * t201 + (t161 * t260 - t28 * t318 + t56 * t67) * t200 + (t209 * t28 + t67 * t86) * t199 + t209 * t1 + t86 * t4) * t179;];
tauX  = t49;
