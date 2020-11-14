% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G1A0
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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:24
% EndTime: 2020-08-06 19:01:31
% DurationCPUTime: 7.54s
% Computational Cost: add. (22827->474), mult. (36156->787), div. (4953->6), fcn. (32256->18), ass. (0->327)
t224 = pkin(1) + pkin(2);
t207 = sin(qJ(2,3));
t221 = xDP(3);
t222 = xDP(2);
t223 = xDP(1);
t226 = 0.1e1 / qJ(3,3);
t213 = cos(qJ(2,3));
t310 = t213 * t224;
t333 = qJ(3,3) * t207;
t155 = t310 + t333;
t149 = 0.1e1 / t155;
t325 = t149 * t213;
t201 = legFrame(3,3);
t167 = sin(t201);
t170 = cos(t201);
t208 = sin(qJ(1,3));
t214 = cos(qJ(1,3));
t118 = -t167 * t208 + t170 * t214;
t179 = t214 * pkin(4);
t143 = t208 * t333 + t179;
t364 = pkin(4) * t208;
t144 = t214 * t333 - t364;
t76 = t118 * t310 - t143 * t167 + t144 * t170;
t119 = t167 * t214 + t170 * t208;
t77 = t119 * t310 + t143 * t170 + t144 * t167;
t58 = (t207 * t221 + (t222 * t77 + t223 * t76) * t325) * t226;
t346 = t224 * t58;
t317 = t207 * t224;
t332 = qJ(3,3) * t213;
t249 = -t317 + t332;
t112 = t155 * t208 + t179;
t115 = t155 * t214 - t364;
t94 = t112 * t170 + t115 * t167;
t342 = t226 * t94;
t91 = -t112 * t167 + t115 * t170;
t343 = t226 * t91;
t61 = -t221 * t226 * t249 + t222 * t342 + t223 * t343;
t382 = t61 - t346;
t209 = sin(qJ(2,2));
t228 = 0.1e1 / qJ(3,2);
t215 = cos(qJ(2,2));
t308 = t215 * t224;
t335 = qJ(3,2) * t209;
t156 = t308 + t335;
t150 = 0.1e1 / t156;
t323 = t150 * t215;
t202 = legFrame(2,3);
t168 = sin(t202);
t171 = cos(t202);
t210 = sin(qJ(1,2));
t216 = cos(qJ(1,2));
t120 = -t168 * t210 + t171 * t216;
t180 = t216 * pkin(4);
t145 = t210 * t335 + t180;
t362 = pkin(4) * t210;
t146 = t216 * t335 - t362;
t78 = t120 * t308 - t145 * t168 + t146 * t171;
t121 = t168 * t216 + t171 * t210;
t79 = t121 * t308 + t145 * t171 + t146 * t168;
t59 = (t209 * t221 + (t222 * t79 + t223 * t78) * t323) * t228;
t345 = t224 * t59;
t314 = t209 * t224;
t334 = qJ(3,2) * t215;
t250 = -t314 + t334;
t113 = t156 * t210 + t180;
t116 = t156 * t216 - t362;
t95 = t113 * t171 + t116 * t168;
t340 = t228 * t95;
t92 = -t113 * t168 + t116 * t171;
t341 = t228 * t92;
t62 = -t221 * t228 * t250 + t222 * t340 + t223 * t341;
t381 = t62 - t345;
t211 = sin(qJ(2,1));
t230 = 0.1e1 / qJ(3,1);
t217 = cos(qJ(2,1));
t306 = t217 * t224;
t337 = qJ(3,1) * t211;
t157 = t306 + t337;
t151 = 0.1e1 / t157;
t321 = t151 * t217;
t203 = legFrame(1,3);
t169 = sin(t203);
t172 = cos(t203);
t212 = sin(qJ(1,1));
t218 = cos(qJ(1,1));
t122 = -t169 * t212 + t172 * t218;
t181 = t218 * pkin(4);
t147 = t212 * t337 + t181;
t360 = pkin(4) * t212;
t148 = t218 * t337 - t360;
t80 = t122 * t306 - t147 * t169 + t148 * t172;
t123 = t169 * t218 + t172 * t212;
t81 = t123 * t306 + t147 * t172 + t148 * t169;
t60 = (t211 * t221 + (t222 * t81 + t223 * t80) * t321) * t230;
t344 = t224 * t60;
t311 = t211 * t224;
t336 = qJ(3,1) * t217;
t251 = -t311 + t336;
t114 = t157 * t212 + t181;
t117 = t157 * t218 - t360;
t96 = t114 * t172 + t117 * t169;
t338 = t230 * t96;
t93 = -t114 * t169 + t117 * t172;
t339 = t230 * t93;
t63 = -t221 * t230 * t251 + t222 * t338 + t223 * t339;
t380 = t63 - t344;
t225 = qJ(3,3) ^ 2;
t376 = 2 * rSges(3,3);
t379 = qJ(3,3) * t376 + t225;
t227 = qJ(3,2) ^ 2;
t378 = qJ(3,2) * t376 + t227;
t229 = qJ(3,1) ^ 2;
t377 = qJ(3,1) * t376 + t229;
t375 = m(1) * rSges(1,1);
t220 = m(2) * rSges(2,2);
t374 = m(2) * rSges(2,3);
t373 = (rSges(3,2) * m(3));
t372 = pkin(4) * t58;
t371 = pkin(4) * t59;
t370 = pkin(4) * t60;
t219 = pkin(1) + rSges(3,1);
t198 = -qJ(3,3) - rSges(3,3);
t369 = m(3) * t198;
t199 = -qJ(3,2) - rSges(3,3);
t368 = m(3) * t199;
t200 = -qJ(3,1) - rSges(3,3);
t367 = m(3) * t200;
t366 = m(3) * t219;
t365 = pkin(4) * t207;
t363 = pkin(4) * t209;
t361 = pkin(4) * t211;
t165 = m(2) * rSges(2,1) + t366;
t359 = g(3) * t165;
t358 = rSges(3,2) * t207;
t357 = rSges(3,2) * t209;
t356 = rSges(3,2) * t211;
t194 = t213 ^ 2;
t73 = (t118 * t222 - t119 * t223) * t149;
t355 = t194 * t73;
t195 = t215 ^ 2;
t74 = (t120 * t222 - t121 * t223) * t150;
t354 = t195 * t74;
t196 = t217 ^ 2;
t75 = (t122 * t222 - t123 * t223) * t151;
t353 = t196 * t75;
t352 = t198 * t61;
t351 = t199 * t62;
t350 = t200 * t63;
t349 = t213 * t73;
t348 = t215 * t74;
t347 = t217 * t75;
t331 = t118 * t149;
t330 = t119 * t149;
t329 = t120 * t150;
t328 = t121 * t150;
t327 = t122 * t151;
t326 = t123 * t151;
t319 = t207 * t213;
t318 = t207 * t219;
t316 = t209 * t215;
t315 = t209 * t219;
t313 = t211 * t217;
t312 = t211 * t219;
t309 = t213 * t226;
t307 = t215 * t228;
t305 = t217 * t230;
t304 = t219 * t226;
t303 = t219 * t228;
t302 = t219 * t230;
t301 = -Icges(2,1) - Icges(3,1);
t231 = rSges(3,3) ^ 2;
t300 = rSges(3,2) ^ 2 + t231;
t299 = -2 * t373;
t298 = -0.2e1 * qJ(3,1) * (t75 * t311 - t370 / 0.2e1);
t297 = -0.2e1 * qJ(3,2) * (t74 * t314 - t371 / 0.2e1);
t296 = -0.2e1 * qJ(3,3) * (t73 * t317 - t372 / 0.2e1);
t295 = rSges(3,3) + t219;
t294 = rSges(3,3) - t219;
t293 = m(3) * t358;
t292 = m(3) * t357;
t291 = m(3) * t356;
t290 = m(3) * t352;
t289 = m(3) * t351;
t288 = m(3) * t350;
t287 = m(3) * t304;
t286 = m(3) * t303;
t285 = m(3) * t302;
t284 = t76 * t325;
t283 = t77 * t325;
t233 = rSges(2,2) ^ 2;
t235 = rSges(2,1) ^ 2;
t256 = (t233 + t235) * m(2) + Icges(3,2) + Icges(2,3);
t238 = (pkin(1) ^ 2);
t257 = t231 + t238 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t103 = (t257 + t379) * m(3) + t256;
t97 = (t103 * t207 + t249 * t366) * t226;
t282 = t97 * t325;
t281 = t78 * t323;
t280 = t79 * t323;
t104 = (t257 + t378) * m(3) + t256;
t98 = (t104 * t209 + t250 * t366) * t228;
t279 = t98 * t323;
t278 = t80 * t321;
t277 = t81 * t321;
t105 = (t257 + t377) * m(3) + t256;
t99 = (t105 * t211 + t251 * t366) * t230;
t276 = t99 * t321;
t262 = -rSges(2,2) * t374 + Icges(2,6) - Icges(3,6);
t139 = -rSges(3,2) * t369 + t262;
t142 = rSges(2,1) * t374 + rSges(3,2) * t366 - Icges(3,4) - Icges(2,5);
t100 = t139 * t213 - t207 * t142;
t275 = t100 * t309;
t140 = -rSges(3,2) * t368 + t262;
t101 = t140 * t215 - t209 * t142;
t274 = t101 * t307;
t141 = -rSges(3,2) * t367 + t262;
t102 = t141 * t217 - t211 * t142;
t273 = t102 * t305;
t272 = t103 * t309;
t271 = t104 * t307;
t270 = t105 * t305;
t269 = (qJ(3,3) + t224) * (-qJ(3,3) + t224) * t194;
t268 = (qJ(3,2) + t224) * (-qJ(3,2) + t224) * t195;
t267 = (qJ(3,1) + t224) * (-qJ(3,1) + t224) * t196;
t266 = t213 * t304;
t265 = t215 * t303;
t264 = t217 * t302;
t263 = -rSges(2,1) * t220 + Icges(2,4) - Icges(3,5);
t261 = -t238 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t159 = t220 + t369;
t160 = t220 + t368;
t161 = t220 + t367;
t260 = t226 * t293;
t259 = t228 * t292;
t258 = t230 * t291;
t252 = Icges(2,2) + Icges(3,3) + (-t233 + t235) * m(2) + t301;
t133 = -g(1) * t167 + g(2) * t170;
t136 = g(1) * t170 + g(2) * t167;
t248 = t133 * t208 + t136 * t214;
t134 = -g(1) * t168 + g(2) * t171;
t137 = g(1) * t171 + g(2) * t168;
t247 = t134 * t210 + t137 * t216;
t135 = -g(1) * t169 + g(2) * t172;
t138 = g(1) * t172 + g(2) * t169;
t246 = t135 * t212 + t138 * t218;
t245 = -t159 * t207 + t165 * t213 + t375;
t244 = -t160 * t209 + t165 * t215 + t375;
t243 = -t161 * t211 + t165 * t217 + t375;
t242 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,3) ^ 2 + t233) * m(2) + Icges(1,3) - t301;
t236 = pkin(4) ^ 2;
t206 = xDDP(1);
t205 = xDDP(2);
t204 = xDDP(3);
t175 = t229 + t236;
t174 = t227 + t236;
t173 = t225 + t236;
t158 = m(1) * rSges(1,2) - t373 - t374;
t132 = -t200 * t366 + t263;
t131 = -t199 * t366 + t263;
t130 = -t198 * t366 + t263;
t111 = (-t251 - t312) * t230 * m(3);
t110 = (-t250 - t315) * t228 * m(3);
t109 = (-t249 - t318) * t226 * m(3);
t108 = -(qJ(3,1) + t295) * (qJ(3,1) + t294) * m(3) + t252;
t107 = -(qJ(3,2) + t295) * (qJ(3,2) + t294) * m(3) + t252;
t106 = -(qJ(3,3) + t295) * (qJ(3,3) + t294) * m(3) + t252;
t84 = (-t251 * t373 + t102) * t230 * t211;
t83 = (-t250 * t373 + t101) * t228 * t209;
t82 = (-t249 * t373 + t100) * t226 * t207;
t72 = t75 ^ 2;
t71 = t74 ^ 2;
t70 = t73 ^ 2;
t69 = t75 * t361;
t68 = t74 * t363;
t67 = t73 * t365;
t66 = t108 * t196 + 0.2e1 * t132 * t313 + (t300 + t377) * m(3) + t242;
t65 = t107 * t195 + 0.2e1 * t131 * t316 + (t300 + t378) * m(3) + t242;
t64 = t106 * t194 + 0.2e1 * t130 * t319 + (t300 + t379) * m(3) + t242;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t48 = (t339 + (-t123 * t356 - t80 * t264) * t151) * m(3);
t47 = (t338 + (t122 * t356 - t81 * t264) * t151) * m(3);
t46 = (t341 + (-t121 * t357 - t78 * t265) * t150) * m(3);
t45 = (t340 + (t120 * t357 - t79 * t265) * t150) * m(3);
t44 = (t343 + (-t119 * t358 - t76 * t266) * t149) * m(3);
t43 = (t342 + (t118 * t358 - t77 * t266) * t149) * m(3);
t42 = t132 * t60;
t41 = t131 * t59;
t40 = t130 * t58;
t39 = -t93 * t285 + (-t102 * t123 + t80 * t270) * t151;
t38 = -t96 * t285 + (t102 * t122 + t81 * t270) * t151;
t37 = -t92 * t286 + (-t101 * t121 + t78 * t271) * t150;
t36 = -t95 * t286 + (t101 * t120 + t79 * t271) * t150;
t35 = -t91 * t287 + (-t100 * t119 + t76 * t272) * t149;
t34 = -t94 * t287 + (t100 * t118 + t77 * t272) * t149;
t33 = t69 + t344;
t32 = t68 + t345;
t31 = t67 + t346;
t27 = t93 * t258 + (-t123 * t66 + t80 * t273) * t151;
t26 = t96 * t258 + (t122 * t66 + t81 * t273) * t151;
t25 = t92 * t259 + (-t121 * t65 + t78 * t274) * t150;
t24 = t95 * t259 + (t120 * t65 + t79 * t274) * t150;
t23 = t91 * t260 + (-t119 * t64 + t76 * t275) * t149;
t22 = t94 * t260 + (t118 * t64 + t77 * t275) * t149;
t21 = t380 * t361;
t20 = t381 * t363;
t19 = t382 * t365;
t18 = (-pkin(4) * t75 + (t251 + t336) * t60 + (0.2e1 * t63 - t344) * t211) * t75 * t151;
t17 = (-pkin(4) * t74 + (t250 + t334) * t59 + (0.2e1 * t62 - t345) * t209) * t74 * t150;
t16 = (-pkin(4) * t73 + (t249 + t332) * t58 + (0.2e1 * t61 - t346) * t207) * t73 * t149;
t15 = (((-t229 + t261) * t60 + t224 * t63) * t60 + t33 * t63 + (t217 * t298 + t21 + (-t175 - t267) * t75 + t251 * t370) * t75) * t230;
t14 = (((-t227 + t261) * t59 + t224 * t62) * t59 + t32 * t62 + (t215 * t297 + t20 + (-t174 - t268) * t74 + t250 * t371) * t74) * t228;
t13 = (((-t225 + t261) * t58 + t224 * t61) * t58 + t31 * t61 + (t213 * t296 + t19 + (-t173 - t269) * t73 + t249 * t372) * t73) * t226;
t12 = ((-t267 * t347 + t196 * t298 + (-t175 * t75 + t21) * t217) * t75 + (-(t69 - t380) * t306 + (pkin(4) * t353 + t211 * t380) * qJ(3,1)) * t60 + (t217 * t33 + t60 * t337) * t63) * t151 * t230;
t11 = ((-t268 * t348 + t195 * t297 + (-t174 * t74 + t20) * t215) * t74 + (-(t68 - t381) * t308 + (pkin(4) * t354 + t209 * t381) * qJ(3,2)) * t59 + (t215 * t32 + t59 * t335) * t62) * t150 * t228;
t10 = ((-t269 * t349 + t194 * t296 + (-t173 * t73 + t19) * t213) * t73 + (-(t67 - t382) * t310 + (pkin(4) * t355 + t207 * t382) * qJ(3,3)) * t58 + (t213 * t31 + t58 * t333) * t61) * t149 * t226;
t9 = t12 * t366 - t18 * t291 + (-t15 - (t196 * t72 - t57 - t72) * t200 - t211 * t246 + (-t72 * t312 + g(3)) * t217) * m(3);
t8 = t11 * t366 - t17 * t292 + (-t14 - (t195 * t71 - t56 - t71) * t199 - t209 * t247 + (-t71 * t315 + g(3)) * t215) * m(3);
t7 = t10 * t366 - t16 * t293 + (-t13 - (t194 * t70 - t55 - t70) * t198 - t207 * t248 + (-t70 * t318 + g(3)) * t213) * m(3);
t6 = -t102 * t18 - t105 * t12 + (t246 * t161 - t359) * t217 + t211 * (g(3) * t161 + t246 * t165) + (t15 * t219 - 0.2e1 * t60 * t350) * m(3) + (t108 * t313 - 0.2e1 * t132 * t196 + t132) * t72;
t5 = -t101 * t17 - t104 * t11 + (t247 * t160 - t359) * t215 + t209 * (g(3) * t160 + t247 * t165) + (t14 * t219 - 0.2e1 * t59 * t351) * m(3) + (t107 * t316 - 0.2e1 * t131 * t195 + t131) * t71;
t4 = -t100 * t16 - t103 * t10 + (t248 * t159 - t359) * t213 + t207 * (g(3) * t159 + t248 * t165) + (t13 * t219 - 0.2e1 * t58 * t352) * m(3) + (t106 * t319 - 0.2e1 * t130 * t194 + t130) * t70;
t3 = -t66 * t18 - t102 * t12 - 0.4e1 * (-t42 - t288 / 0.2e1) * t353 - t60 * (t142 * t60 + t63 * t299) * t217 + 0.2e1 * (-t42 - t288) * t75 + (-t15 * t373 - 0.2e1 * (t108 * t60 - t63 * t366) * t347 - t57 * t141) * t211 + (t158 * t218 + t243 * t212) * t138 + (t158 * t212 - t243 * t218) * t135;
t2 = -t65 * t17 - t101 * t11 - 0.4e1 * (-t41 - t289 / 0.2e1) * t354 - t59 * (t142 * t59 + t62 * t299) * t215 + 0.2e1 * (-t41 - t289) * t74 + (-t14 * t373 - 0.2e1 * (t107 * t59 - t62 * t366) * t348 - t56 * t140) * t209 + (t158 * t216 + t244 * t210) * t137 + (t158 * t210 - t244 * t216) * t134;
t1 = -t64 * t16 - t100 * t10 - 0.4e1 * (-t40 - t290 / 0.2e1) * t355 - t58 * (t142 * t58 + t61 * t299) * t213 + 0.2e1 * (-t40 - t290) * t73 + (-t13 * t373 - 0.2e1 * (t106 * t58 - t61 * t366) * t349 - t55 * t139) * t207 + (t158 * t214 + t245 * t208) * t136 + (t158 * t208 - t245 * t214) * t133;
t28 = [-t1 * t330 - t2 * t328 - t3 * t326 - m(4) * g(1) + (t23 * t331 + t25 * t329 + t27 * t327) * t205 + (-t23 * t330 - t25 * t328 - t27 * t326 + m(4)) * t206 + ((t39 * t278 + t48 * t93) * t206 + (t39 * t277 + t48 * t96) * t205 + (t211 * t39 - t251 * t48) * t204 + t6 * t278 + t93 * t9) * t230 + ((t37 * t281 + t46 * t92) * t206 + (t37 * t280 + t46 * t95) * t205 + (t209 * t37 - t250 * t46) * t204 + t5 * t281 + t92 * t8) * t228 + ((t35 * t284 + t44 * t91) * t206 + (t35 * t283 + t44 * t94) * t205 + (t207 * t35 - t249 * t44) * t204 + t4 * t284 + t91 * t7) * t226; t1 * t331 + t2 * t329 + t3 * t327 - m(4) * g(2) + (-t22 * t330 - t24 * t328 - t26 * t326) * t206 + (t22 * t331 + t24 * t329 + t26 * t327 + m(4)) * t205 + ((t278 * t38 + t47 * t93) * t206 + (t277 * t38 + t47 * t96) * t205 + (t211 * t38 - t251 * t47) * t204 + t6 * t277 + t96 * t9) * t230 + ((t281 * t36 + t45 * t92) * t206 + (t280 * t36 + t45 * t95) * t205 + (t209 * t36 - t250 * t45) * t204 + t5 * t280 + t95 * t8) * t228 + ((t284 * t34 + t43 * t91) * t206 + (t283 * t34 + t43 * t94) * t205 + (t207 * t34 - t249 * t43) * t204 + t4 * t283 + t94 * t7) * t226; (-g(3) + t204) * m(4) + (-t84 * t326 - t83 * t328 - t82 * t330) * t206 + (t84 * t327 + t83 * t329 + t82 * t331) * t205 + ((t111 * t93 + t276 * t80) * t206 + (t111 * t96 + t276 * t81) * t205 + (-t111 * t251 + t211 * t99) * t204 + t211 * t6 - t251 * t9) * t230 + ((t110 * t92 + t279 * t78) * t206 + (t110 * t95 + t279 * t79) * t205 + (-t110 * t250 + t209 * t98) * t204 + t209 * t5 - t250 * t8) * t228 + ((t109 * t91 + t282 * t76) * t206 + (t109 * t94 + t282 * t77) * t205 + (-t109 * t249 + t207 * t97) * t204 + t207 * t4 - t249 * t7) * t226;];
tauX  = t28;
