% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G3A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:09:40
% EndTime: 2020-08-06 19:09:47
% DurationCPUTime: 7.46s
% Computational Cost: add. (31194->495), mult. (46491->814), div. (6867->6), fcn. (37389->18), ass. (0->321)
t207 = pkin(1) + pkin(2);
t196 = cos(qJ(2,3));
t275 = t196 * t207;
t190 = sin(qJ(2,3));
t313 = qJ(3,3) * t190;
t146 = t275 + t313;
t191 = sin(qJ(1,3));
t197 = cos(qJ(1,3));
t115 = pkin(4) * t197 + t146 * t191;
t204 = xDP(3);
t205 = xDP(2);
t206 = xDP(1);
t140 = 0.1e1 / t146;
t209 = 0.1e1 / qJ(3,3);
t299 = t140 * t209;
t283 = t190 * t197;
t350 = pkin(4) * t191;
t137 = qJ(3,3) * t283 - t350;
t184 = legFrame(3,2);
t161 = sin(t184);
t164 = cos(t184);
t177 = t196 ^ 2;
t273 = t197 * t207;
t282 = t190 * t207;
t296 = t161 * qJ(3,3);
t85 = (t164 * t273 - t296) * t177 + (t137 * t164 + t161 * t282) * t196 + t296;
t290 = t164 * qJ(3,3);
t88 = (-t161 * t273 - t290) * t177 + (-t137 * t161 + t164 * t282) * t196 + t290;
t52 = (-t115 * t196 * t204 + t205 * t88 + t206 * t85) * t299;
t332 = t207 * t52;
t198 = cos(qJ(2,2));
t272 = t198 * t207;
t192 = sin(qJ(2,2));
t315 = qJ(3,2) * t192;
t147 = t272 + t315;
t193 = sin(qJ(1,2));
t199 = cos(qJ(1,2));
t116 = pkin(4) * t199 + t147 * t193;
t141 = 0.1e1 / t147;
t211 = 0.1e1 / qJ(3,2);
t298 = t141 * t211;
t280 = t192 * t199;
t349 = pkin(4) * t193;
t138 = qJ(3,2) * t280 - t349;
t185 = legFrame(2,2);
t162 = sin(t185);
t165 = cos(t185);
t178 = t198 ^ 2;
t270 = t199 * t207;
t279 = t192 * t207;
t294 = t162 * qJ(3,2);
t86 = (t165 * t270 - t294) * t178 + (t138 * t165 + t162 * t279) * t198 + t294;
t288 = t165 * qJ(3,2);
t89 = (-t162 * t270 - t288) * t178 + (-t138 * t162 + t165 * t279) * t198 + t288;
t53 = (-t116 * t198 * t204 + t205 * t89 + t206 * t86) * t298;
t331 = t207 * t53;
t200 = cos(qJ(2,1));
t269 = t200 * t207;
t194 = sin(qJ(2,1));
t317 = qJ(3,1) * t194;
t148 = t269 + t317;
t195 = sin(qJ(1,1));
t201 = cos(qJ(1,1));
t117 = pkin(4) * t201 + t148 * t195;
t142 = 0.1e1 / t148;
t213 = 0.1e1 / qJ(3,1);
t297 = t142 * t213;
t277 = t194 * t201;
t348 = pkin(4) * t195;
t139 = qJ(3,1) * t277 - t348;
t186 = legFrame(1,2);
t163 = sin(t186);
t166 = cos(t186);
t179 = t200 ^ 2;
t267 = t201 * t207;
t276 = t194 * t207;
t292 = t163 * qJ(3,1);
t87 = (t166 * t267 - t292) * t179 + (t139 * t166 + t163 * t276) * t200 + t292;
t286 = t166 * qJ(3,1);
t90 = (-t163 * t267 - t286) * t179 + (-t139 * t163 + t166 * t276) * t200 + t286;
t54 = (-t117 * t200 * t204 + t205 * t90 + t206 * t87) * t297;
t330 = t207 * t54;
t367 = t146 * t197 - t350;
t366 = t147 * t199 - t349;
t365 = t148 * t201 - t348;
t208 = qJ(3,3) ^ 2;
t361 = 2 * rSges(3,3);
t364 = qJ(3,3) * t361 + t208;
t210 = qJ(3,2) ^ 2;
t363 = qJ(3,2) * t361 + t210;
t212 = qJ(3,1) ^ 2;
t362 = qJ(3,1) * t361 + t212;
t360 = m(1) * rSges(1,1);
t203 = m(2) * rSges(2,2);
t359 = m(2) * rSges(2,3);
t358 = (rSges(3,2) * m(3));
t82 = (-t197 * t204 + (t161 * t205 - t164 * t206) * t191) * t140;
t357 = pkin(4) * t82;
t83 = (-t199 * t204 + (t162 * t205 - t165 * t206) * t193) * t141;
t356 = pkin(4) * t83;
t84 = (-t201 * t204 + (t163 * t205 - t166 * t206) * t195) * t142;
t355 = pkin(4) * t84;
t202 = pkin(1) + rSges(3,1);
t181 = -qJ(3,3) - rSges(3,3);
t354 = m(3) * t181;
t182 = -qJ(3,2) - rSges(3,3);
t353 = m(3) * t182;
t183 = -qJ(3,1) - rSges(3,3);
t352 = m(3) * t183;
t351 = m(3) * t202;
t344 = rSges(3,2) * t190;
t343 = rSges(3,2) * t192;
t342 = rSges(3,2) * t194;
t341 = t140 * t82;
t340 = t141 * t83;
t339 = t142 * t84;
t302 = t115 * t209;
t312 = qJ(3,3) * t196;
t226 = -t282 + t312;
t100 = -t367 * t161 - t226 * t164;
t311 = t100 * t209;
t97 = -t226 * t161 + t367 * t164;
t327 = t209 * t97;
t67 = -t204 * t302 + t205 * t311 + t206 * t327;
t338 = t181 * t67;
t301 = t116 * t211;
t314 = qJ(3,2) * t198;
t227 = -t279 + t314;
t101 = -t366 * t162 - t227 * t165;
t310 = t101 * t211;
t98 = -t227 * t162 + t366 * t165;
t324 = t211 * t98;
t68 = -t204 * t301 + t205 * t310 + t206 * t324;
t337 = t182 * t68;
t300 = t117 * t213;
t316 = qJ(3,1) * t200;
t228 = -t276 + t316;
t102 = -t365 * t163 - t228 * t166;
t309 = t102 * t213;
t99 = -t228 * t163 + t365 * t166;
t321 = t213 * t99;
t69 = -t204 * t300 + t205 * t309 + t206 * t321;
t336 = t183 * t69;
t22 = -t67 + t332;
t335 = t190 * t22;
t23 = -t68 + t331;
t334 = t192 * t23;
t24 = -t69 + t330;
t333 = t194 * t24;
t329 = t209 * t85;
t328 = t209 * t88;
t326 = t211 * t86;
t325 = t211 * t89;
t323 = t213 * t87;
t322 = t213 * t90;
t320 = t82 * t177;
t319 = t83 * t178;
t318 = t84 * t179;
t239 = -rSges(2,2) * t359 + Icges(2,6) - Icges(3,6);
t127 = -rSges(3,2) * t354 + t239;
t136 = rSges(2,1) * t359 + rSges(3,2) * t351 - Icges(3,4) - Icges(2,5);
t103 = t127 * t196 - t136 * t190;
t308 = t103 * t209;
t128 = -rSges(3,2) * t353 + t239;
t104 = t128 * t198 - t136 * t192;
t307 = t104 * t211;
t129 = -rSges(3,2) * t352 + t239;
t105 = t129 * t200 - t136 * t194;
t306 = t105 * t213;
t216 = rSges(2,2) ^ 2;
t218 = rSges(2,1) ^ 2;
t233 = (t216 + t218) * m(2) + Icges(3,2) + Icges(2,3);
t214 = rSges(3,3) ^ 2;
t221 = (pkin(1) ^ 2);
t234 = t214 + t221 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t106 = (t234 + t364) * m(3) + t233;
t305 = t106 * t209;
t107 = (t234 + t363) * m(3) + t233;
t304 = t107 * t211;
t108 = (t234 + t362) * m(3) + t233;
t303 = t108 * t213;
t295 = t161 * t191;
t293 = t162 * t193;
t291 = t163 * t195;
t289 = t164 * t191;
t287 = t165 * t193;
t285 = t166 * t195;
t284 = t190 * t196;
t281 = t192 * t198;
t278 = t194 * t200;
t274 = t196 * t209;
t271 = t198 * t211;
t268 = t200 * t213;
t266 = t202 * t209;
t265 = t202 * t211;
t264 = t202 * t213;
t263 = -Icges(2,1) - Icges(3,1);
t262 = rSges(3,2) ^ 2 + t214;
t261 = -2 * t358;
t260 = rSges(3,3) - t202;
t259 = rSges(3,3) + t202;
t258 = m(3) * t344;
t257 = m(3) * t343;
t256 = m(3) * t342;
t255 = m(3) * t338;
t254 = m(3) * t337;
t253 = m(3) * t336;
t252 = 0.2e1 * t284;
t251 = 0.2e1 * t281;
t250 = 0.2e1 * t278;
t249 = m(3) * t266;
t248 = m(3) * t265;
t247 = m(3) * t264;
t246 = t191 * t344;
t245 = t193 * t343;
t244 = t195 * t342;
t243 = t115 * t274;
t242 = t116 * t271;
t241 = t117 * t268;
t240 = -rSges(2,1) * t203 + Icges(2,4) - Icges(3,5);
t238 = -t221 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t153 = t203 + t354;
t154 = t203 + t353;
t155 = t203 + t352;
t237 = t209 * t258;
t236 = t211 * t257;
t235 = t213 * t256;
t133 = g(1) * t164 - g(2) * t161;
t232 = -g(3) * t191 + t133 * t197;
t134 = g(1) * t165 - g(2) * t162;
t231 = -g(3) * t193 + t134 * t199;
t135 = g(1) * t166 - g(2) * t163;
t230 = -g(3) * t195 + t135 * t201;
t229 = Icges(2,2) + Icges(3,3) + (-t216 + t218) * m(2) + t263;
t159 = m(2) * rSges(2,1) + t351;
t225 = -t153 * t190 + t159 * t196;
t224 = -t154 * t192 + t159 * t198;
t223 = -t155 * t194 + t159 * t200;
t222 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,3) ^ 2 + t216) * m(2) + Icges(1,3) - t263;
t219 = pkin(4) ^ 2;
t189 = xDDP(1);
t188 = xDDP(2);
t187 = xDDP(3);
t175 = g(3) * t360;
t152 = m(1) * rSges(1,2) - t358 - t359;
t132 = g(1) * t163 + g(2) * t166;
t131 = g(1) * t162 + g(2) * t165;
t130 = g(1) * t161 + g(2) * t164;
t126 = -t183 * t351 + t240;
t125 = -t182 * t351 + t240;
t124 = -t181 * t351 + t240;
t114 = -(qJ(3,1) + t259) * (qJ(3,1) + t260) * m(3) + t229;
t113 = -(qJ(3,2) + t259) * (qJ(3,2) + t260) * m(3) + t229;
t112 = -(qJ(3,3) + t259) * (qJ(3,3) + t260) * m(3) + t229;
t81 = t84 ^ 2;
t80 = t83 ^ 2;
t79 = t82 ^ 2;
t78 = t194 * t355;
t77 = t192 * t356;
t76 = t190 * t357;
t75 = t114 * t179 + t126 * t250 + (t262 + t362) * m(3) + t222;
t74 = t113 * t178 + t125 * t251 + (t262 + t363) * m(3) + t222;
t73 = t112 * t177 + t124 * t252 + (t262 + t364) * m(3) + t222;
t72 = (-t300 + (-rSges(3,2) * t277 + t202 * t241) * t142) * m(3);
t71 = (-t301 + (-rSges(3,2) * t280 + t202 * t242) * t141) * m(3);
t70 = (-t302 + (-rSges(3,2) * t283 + t202 * t243) * t140) * m(3);
t66 = t117 * t247 + (-t105 * t201 - t108 * t241) * t142;
t65 = t116 * t248 + (-t104 * t199 - t107 * t242) * t141;
t64 = t115 * t249 + (-t103 * t197 - t106 * t243) * t140;
t63 = (t321 + (-t166 * t244 - t87 * t264) * t142) * m(3);
t62 = (t324 + (-t165 * t245 - t86 * t265) * t141) * m(3);
t61 = (t327 + (-t164 * t246 - t85 * t266) * t140) * m(3);
t60 = (t309 + (t163 * t244 - t90 * t264) * t142) * m(3);
t59 = (t310 + (t162 * t245 - t89 * t265) * t141) * m(3);
t58 = (t311 + (t161 * t246 - t88 * t266) * t140) * m(3);
t57 = -t117 * t235 + (-t105 * t241 - t201 * t75) * t142;
t56 = -t116 * t236 + (-t104 * t242 - t199 * t74) * t141;
t55 = -t115 * t237 + (-t103 * t243 - t197 * t73) * t140;
t51 = t54 ^ 2;
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t42 = -t99 * t247 + (-t105 * t285 + t87 * t303) * t142;
t41 = -t98 * t248 + (-t104 * t287 + t86 * t304) * t141;
t40 = -t97 * t249 + (-t103 * t289 + t85 * t305) * t140;
t39 = -t102 * t247 + (t105 * t291 + t90 * t303) * t142;
t38 = -t101 * t248 + (t104 * t293 + t89 * t304) * t141;
t37 = -t100 * t249 + (t103 * t295 + t88 * t305) * t140;
t36 = t126 * t54;
t35 = t125 * t53;
t34 = t124 * t52;
t33 = t99 * t235 + (-t75 * t285 + t87 * t306) * t142;
t32 = t98 * t236 + (-t74 * t287 + t86 * t307) * t141;
t31 = t97 * t237 + (-t73 * t289 + t85 * t308) * t140;
t30 = t102 * t235 + (t75 * t291 + t90 * t306) * t142;
t29 = t101 * t236 + (t74 * t293 + t89 * t307) * t141;
t28 = t100 * t237 + (t73 * t295 + t88 * t308) * t140;
t27 = t78 + t330;
t26 = t77 + t331;
t25 = t76 + t332;
t21 = (-t54 * t316 + t333) * pkin(4) + ((qJ(3,1) + t207) * (-qJ(3,1) + t207) * t179 + qJ(3,1) * t207 * t250 + t212 + t219) * t84;
t20 = (-t53 * t314 + t334) * pkin(4) + ((qJ(3,2) + t207) * (-qJ(3,2) + t207) * t178 + t207 * qJ(3,2) * t251 + t210 + t219) * t83;
t19 = (-t52 * t312 + t335) * pkin(4) + ((qJ(3,3) + t207) * (-qJ(3,3) + t207) * t177 + t207 * qJ(3,3) * t252 + t208 + t219) * t82;
t18 = (-t355 + (t228 + t316) * t54 + (0.2e1 * t69 - t330) * t194) * t339;
t17 = (-t356 + (t227 + t314) * t53 + (0.2e1 * t68 - t331) * t192) * t340;
t16 = (-t357 + (t226 + t312) * t52 + (0.2e1 * t67 - t332) * t190) * t341;
t15 = (((-t212 + t238) * t54 + t207 * t69) * t54 + t27 * t69 + (t228 * t54 * pkin(4) - t21) * t84) * t213;
t14 = (((-t210 + t238) * t53 + t207 * t68) * t53 + t26 * t68 + (t227 * t53 * pkin(4) - t20) * t83) * t211;
t13 = (((-t208 + t238) * t52 + t207 * t67) * t52 + t25 * t67 + (t226 * t52 * pkin(4) - t19) * t82) * t209;
t12 = -t21 * t268 * t339 + ((-(t78 + t24) * t269 + (pkin(4) * t318 - t333) * qJ(3,1)) * t54 + (t200 * t27 + t54 * t317) * t69) * t297;
t11 = -t20 * t271 * t340 + ((-(t77 + t23) * t272 + (pkin(4) * t319 - t334) * qJ(3,2)) * t53 + (t198 * t26 + t53 * t315) * t68) * t298;
t10 = -t19 * t274 * t341 + ((-(t76 + t22) * t275 + (pkin(4) * t320 - t335) * qJ(3,3)) * t52 + (t196 * t25 + t52 * t313) * t67) * t299;
t9 = t12 * t351 - t18 * t256 + (-t15 - t81 * t202 * t278 - (t81 * t179 - t51 - t81) * t183 + t132 * t200 - t230 * t194) * m(3);
t8 = t11 * t351 - t17 * t257 + (-t14 - t80 * t202 * t281 - (t178 * t80 - t50 - t80) * t182 + t131 * t198 - t231 * t192) * m(3);
t7 = t10 * t351 - t16 * t258 + (-t13 - t79 * t202 * t284 - (t177 * t79 - t49 - t79) * t181 + t130 * t196 - t232 * t190) * m(3);
t6 = -t105 * t18 - t108 * t12 + (-t132 * t159 + t230 * t155) * t200 + (t132 * t155 + t230 * t159) * t194 + (t202 * t15 - 0.2e1 * t54 * t336) * m(3) + (t114 * t278 - 0.2e1 * t126 * t179 + t126) * t81;
t5 = -t104 * t17 - t107 * t11 + (-t159 * t131 + t231 * t154) * t198 + (t131 * t154 + t231 * t159) * t192 + (t202 * t14 - 0.2e1 * t53 * t337) * m(3) + (t113 * t281 - 0.2e1 * t125 * t178 + t125) * t80;
t4 = -t103 * t16 - t106 * t10 + (-t130 * t159 + t232 * t153) * t196 + (t130 * t153 + t232 * t159) * t190 + (t13 * t202 - 0.2e1 * t52 * t338) * m(3) + (t112 * t284 - 0.2e1 * t124 * t177 + t124) * t79;
t3 = -t75 * t18 - t105 * t12 - 0.4e1 * (-t36 - t253 / 0.2e1) * t318 - t54 * (t136 * t54 + t69 * t261) * t200 + 0.2e1 * (-t36 - t253) * t84 + t175 * t201 + (-t15 * t358 - 0.2e1 * (t114 * t54 - t69 * t351) * t84 * t200 - t51 * t129) * t194 + (-t152 * t195 + t223 * t201) * g(3) + (t152 * t201 + (t223 + t360) * t195) * t135;
t2 = -t74 * t17 - t104 * t11 - 0.4e1 * (-t35 - t254 / 0.2e1) * t319 - t53 * (t136 * t53 + t68 * t261) * t198 + 0.2e1 * (-t35 - t254) * t83 + t175 * t199 + (-t14 * t358 - 0.2e1 * (t113 * t53 - t68 * t351) * t83 * t198 - t50 * t128) * t192 + (-t152 * t193 + t224 * t199) * g(3) + (t152 * t199 + (t224 + t360) * t193) * t134;
t1 = -t73 * t16 - t103 * t10 - 0.4e1 * (-t34 - t255 / 0.2e1) * t320 - t52 * (t136 * t52 + t67 * t261) * t196 + 0.2e1 * (-t34 - t255) * t82 + t175 * t197 + (-t13 * t358 - 0.2e1 * (t112 * t52 - t67 * t351) * t82 * t196 - t49 * t127) * t190 + (-t152 * t191 + t225 * t197) * g(3) + (t152 * t197 + (t225 + t360) * t191) * t133;
t43 = [t7 * t327 + t8 * t324 + t9 * t321 - m(4) * g(1) + (t63 * t309 + t62 * t310 + t61 * t311) * t188 + (-t63 * t300 - t62 * t301 - t61 * t302) * t187 + (t63 * t321 + t62 * t324 + t61 * t327 + m(4)) * t189 + ((-t33 * t285 + t42 * t323) * t189 + (t33 * t291 + t42 * t322) * t188 + (-t201 * t33 - t42 * t241) * t187 - t3 * t285 + t6 * t323) * t142 + ((-t32 * t287 + t41 * t326) * t189 + (t32 * t293 + t41 * t325) * t188 + (-t199 * t32 - t41 * t242) * t187 - t2 * t287 + t5 * t326) * t141 + ((-t31 * t289 + t40 * t329) * t189 + (t31 * t295 + t40 * t328) * t188 + (-t197 * t31 - t40 * t243) * t187 - t1 * t289 + t4 * t329) * t140; t7 * t311 + t8 * t310 + t9 * t309 - m(4) * g(2) + (t60 * t321 + t59 * t324 + t58 * t327) * t189 + (-t60 * t300 - t59 * t301 - t58 * t302) * t187 + (t60 * t309 + t59 * t310 + t58 * t311 + m(4)) * t188 + ((-t30 * t285 + t39 * t323) * t189 + (t30 * t291 + t39 * t322) * t188 + (-t201 * t30 - t39 * t241) * t187 + t3 * t291 + t6 * t322) * t142 + ((-t29 * t287 + t38 * t326) * t189 + (t29 * t293 + t38 * t325) * t188 + (-t199 * t29 - t38 * t242) * t187 + t2 * t293 + t5 * t325) * t141 + ((-t28 * t289 + t37 * t329) * t189 + (t28 * t295 + t37 * t328) * t188 + (-t197 * t28 - t37 * t243) * t187 + t1 * t295 + t4 * t328) * t140; -t7 * t302 - t8 * t301 - t9 * t300 - m(4) * g(3) + (t72 * t321 + t71 * t324 + t70 * t327) * t189 + (t72 * t309 + t71 * t310 + t70 * t311) * t188 + (-t72 * t300 - t71 * t301 - t70 * t302 + m(4)) * t187 + ((-t57 * t285 + t66 * t323) * t189 + (t57 * t291 + t66 * t322) * t188 + (-t201 * t57 - t66 * t241) * t187 - t201 * t3 - t6 * t241) * t142 + ((-t56 * t287 + t65 * t326) * t189 + (t56 * t293 + t65 * t325) * t188 + (-t199 * t56 - t65 * t242) * t187 - t199 * t2 - t5 * t242) * t141 + ((-t55 * t289 + t64 * t329) * t189 + (t55 * t295 + t64 * t328) * t188 + (-t197 * t55 - t64 * t243) * t187 - t197 * t1 - t4 * t243) * t140;];
tauX  = t43;
