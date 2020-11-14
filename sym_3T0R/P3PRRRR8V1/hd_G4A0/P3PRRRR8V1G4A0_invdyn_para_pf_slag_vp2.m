% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:25:03
% EndTime: 2020-08-06 17:25:12
% DurationCPUTime: 9.53s
% Computational Cost: add. (43494->483), mult. (113499->964), div. (7479->8), fcn. (142872->34), ass. (0->360)
t241 = cos(qJ(3,1));
t212 = 0.1e1 / t241;
t242 = cos(qJ(2,1));
t236 = sin(qJ(2,1));
t303 = t236 * t241;
t170 = pkin(2) * t303 - pkin(5) * t242;
t214 = sin(pkin(3));
t235 = sin(qJ(3,1));
t216 = cos(pkin(3));
t381 = pkin(2) * t216;
t146 = t170 * t214 + t235 * t381;
t394 = 0.1e1 / t146;
t344 = t394 * t212;
t239 = cos(qJ(3,2));
t210 = 0.1e1 / t239;
t240 = cos(qJ(2,2));
t234 = sin(qJ(2,2));
t305 = t234 * t239;
t169 = pkin(2) * t305 - pkin(5) * t240;
t233 = sin(qJ(3,2));
t145 = t169 * t214 + t233 * t381;
t395 = 0.1e1 / t145;
t345 = t395 * t210;
t231 = sin(qJ(3,3));
t237 = cos(qJ(3,3));
t284 = t237 * mrSges(3,1) - mrSges(3,2) * t231;
t283 = t239 * mrSges(3,1) - mrSges(3,2) * t233;
t282 = t241 * mrSges(3,1) - mrSges(3,2) * t235;
t222 = legFrame(1,1);
t190 = sin(t222);
t196 = cos(t222);
t227 = legFrame(1,2);
t202 = sin(t227);
t205 = cos(t227);
t135 = g(1) * t202 + (-g(2) * t190 + g(3) * t196) * t205;
t219 = legFrame(1,3);
t187 = sin(t219);
t193 = cos(t219);
t329 = t196 * t202;
t332 = t190 * t202;
t375 = g(1) * t205;
t107 = -t187 * t375 + (-t187 * t332 + t193 * t196) * g(2) + (t187 * t329 + t190 * t193) * g(3);
t108 = t193 * t375 + (t187 * t196 + t193 * t332) * g(2) + (t187 * t190 - t193 * t329) * g(3);
t213 = sin(pkin(6));
t215 = cos(pkin(6));
t90 = t107 * t215 - t108 * t213;
t400 = -t135 * t216 + t90 * t214;
t221 = legFrame(2,1);
t189 = sin(t221);
t195 = cos(t221);
t226 = legFrame(2,2);
t201 = sin(t226);
t204 = cos(t226);
t134 = g(1) * t201 + (-g(2) * t189 + g(3) * t195) * t204;
t218 = legFrame(2,3);
t186 = sin(t218);
t192 = cos(t218);
t330 = t195 * t201;
t333 = t189 * t201;
t376 = g(1) * t204;
t105 = -t186 * t376 + (-t186 * t333 + t192 * t195) * g(2) + (t186 * t330 + t189 * t192) * g(3);
t106 = t192 * t376 + (t186 * t195 + t192 * t333) * g(2) + (t186 * t189 - t192 * t330) * g(3);
t89 = t105 * t215 - t106 * t213;
t399 = -t134 * t216 + t89 * t214;
t220 = legFrame(3,1);
t188 = sin(t220);
t194 = cos(t220);
t225 = legFrame(3,2);
t200 = sin(t225);
t203 = cos(t225);
t133 = g(1) * t200 + (-g(2) * t188 + g(3) * t194) * t203;
t217 = legFrame(3,3);
t185 = sin(t217);
t191 = cos(t217);
t331 = t194 * t200;
t334 = t188 * t200;
t377 = g(1) * t203;
t103 = -t185 * t377 + (-t185 * t334 + t191 * t194) * g(2) + (t185 * t331 + t188 * t191) * g(3);
t104 = t191 * t377 + (t185 * t194 + t191 * t334) * g(2) + (t185 * t188 - t191 * t331) * g(3);
t88 = t103 * t215 - t104 * t213;
t398 = -t133 * t216 + t88 * t214;
t397 = 2 * Ifges(3,4);
t232 = sin(qJ(2,3));
t307 = t232 * t237;
t287 = t214 * t307;
t316 = t216 * t231;
t238 = cos(qJ(2,3));
t321 = t214 * t238;
t396 = 0.1e1 / (-pkin(5) * t321 + (t287 + t316) * pkin(2));
t207 = t237 ^ 2;
t393 = 0.2e1 * t207;
t209 = t239 ^ 2;
t392 = 0.2e1 * t209;
t211 = t241 ^ 2;
t391 = 0.2e1 * t211;
t390 = Ifges(3,5) / 0.2e1;
t389 = -Ifges(3,6) / 0.2e1;
t243 = xDP(3);
t244 = xDP(2);
t248 = 0.1e1 / pkin(2);
t208 = 0.1e1 / t237;
t346 = t396 * t208;
t293 = t248 * t346;
t245 = xDP(1);
t328 = t203 * t245;
t150 = t185 * t215 + t191 * t213;
t153 = -t185 * t213 + t191 * t215;
t109 = t150 * t331 + t153 * t188;
t110 = -t150 * t188 + t153 * t331;
t312 = t216 * t238;
t315 = t216 * t232;
t380 = pkin(2) * t237;
t73 = (-t109 * t232 + t110 * t312) * t380 + pkin(5) * (t109 * t238 + t110 * t315);
t115 = t150 * t194 + t153 * t334;
t118 = t150 * t334 - t153 * t194;
t74 = -(t115 * t312 - t118 * t232) * t380 - pkin(5) * (t115 * t315 + t118 * t238);
t97 = (-t150 * t232 + t153 * t312) * t380 + pkin(5) * (t150 * t238 + t153 * t315);
t49 = (t243 * t73 + t244 * t74 - t97 * t328) * t293;
t388 = pkin(2) * t49;
t292 = t248 * t345;
t327 = t204 * t245;
t151 = t186 * t215 + t192 * t213;
t154 = -t186 * t213 + t192 * t215;
t111 = t151 * t330 + t154 * t189;
t112 = -t151 * t189 + t154 * t330;
t311 = t216 * t240;
t314 = t216 * t234;
t379 = pkin(2) * t239;
t75 = (-t111 * t234 + t112 * t311) * t379 + pkin(5) * (t111 * t240 + t112 * t314);
t116 = t151 * t195 + t154 * t333;
t119 = t151 * t333 - t154 * t195;
t76 = -(t116 * t311 - t119 * t234) * t379 - pkin(5) * (t116 * t314 + t119 * t240);
t98 = (-t151 * t234 + t154 * t311) * t379 + pkin(5) * (t151 * t240 + t154 * t314);
t50 = (t243 * t75 + t244 * t76 - t98 * t327) * t292;
t387 = pkin(2) * t50;
t291 = t248 * t344;
t326 = t205 * t245;
t152 = t187 * t215 + t193 * t213;
t155 = -t187 * t213 + t193 * t215;
t113 = t152 * t329 + t155 * t190;
t114 = -t152 * t190 + t155 * t329;
t310 = t216 * t242;
t313 = t216 * t236;
t378 = pkin(2) * t241;
t77 = (-t113 * t236 + t114 * t310) * t378 + pkin(5) * (t113 * t242 + t114 * t313);
t117 = t152 * t196 + t155 * t332;
t120 = t152 * t332 - t155 * t196;
t78 = -(t117 * t310 - t120 * t236) * t378 - pkin(5) * (t117 * t313 + t120 * t242);
t99 = (-t152 * t236 + t155 * t310) * t378 + pkin(5) * (t152 * t242 + t155 * t313);
t51 = (t243 * t77 + t244 * t78 - t99 * t326) * t291;
t386 = pkin(2) * t51;
t224 = mrSges(2,2) - mrSges(3,3);
t385 = t224 / 0.2e1;
t384 = pkin(2) * t207;
t383 = pkin(2) * t209;
t382 = pkin(2) * t211;
t246 = pkin(5) ^ 2;
t247 = pkin(2) ^ 2;
t356 = t233 * t50;
t301 = pkin(2) * t356;
t157 = t213 * t240 + t215 * t314;
t160 = -t213 * t314 + t215 * t240;
t128 = t157 * t192 + t160 * t186;
t320 = t214 * t239;
t101 = t128 * t233 + t154 * t320;
t259 = t157 * t186 - t160 * t192;
t80 = (t128 * t330 - t259 * t189) * t233 + t112 * t320;
t83 = (-t128 * t333 - t259 * t195) * t233 - t116 * t320;
t56 = (-t101 * t327 + t243 * t80 + t244 * t83) * t345;
t374 = (-pkin(5) * t301 + (t209 * t247 + t246) * t56) * t56;
t354 = t235 * t51;
t300 = pkin(2) * t354;
t158 = t213 * t242 + t215 * t313;
t161 = -t213 * t313 + t215 * t242;
t129 = t158 * t193 + t161 * t187;
t318 = t214 * t241;
t102 = t129 * t235 + t155 * t318;
t258 = t158 * t187 - t161 * t193;
t81 = (t129 * t329 - t258 * t190) * t235 + t114 * t318;
t84 = (-t129 * t332 - t258 * t196) * t235 - t117 * t318;
t57 = (-t102 * t326 + t243 * t81 + t244 * t84) * t344;
t373 = (-pkin(5) * t300 + (t211 * t247 + t246) * t57) * t57;
t372 = Ifges(3,1) + Ifges(2,3);
t156 = t213 * t238 + t215 * t315;
t159 = -t213 * t315 + t215 * t238;
t127 = t156 * t191 + t159 * t185;
t322 = t214 * t237;
t100 = t127 * t231 + t153 * t322;
t168 = pkin(2) * t307 - pkin(5) * t238;
t144 = pkin(2) * t316 + t168 * t214;
t141 = 0.1e1 / t144;
t260 = t156 * t185 - t159 * t191;
t79 = (t127 * t331 - t260 * t188) * t231 + t110 * t322;
t82 = (-t127 * t334 - t260 * t194) * t231 - t115 * t322;
t55 = (-t100 * t328 + t243 * t79 + t244 * t82) * t208 * t141;
t368 = t396 * t55;
t367 = t203 * t97;
t366 = t204 * t98;
t365 = t205 * t99;
t364 = t208 * t79;
t363 = t208 * t82;
t362 = t210 * t80;
t361 = t210 * t83;
t360 = t212 * t81;
t359 = t212 * t84;
t358 = t231 * t49;
t357 = t231 * t55;
t355 = t233 * t56;
t353 = t235 * t57;
t352 = t100 * t203;
t351 = t101 * t204;
t350 = t102 * t205;
t174 = mrSges(2,1) + t284;
t343 = (t174 * t238 - t224 * t232) * t214;
t175 = mrSges(2,1) + t283;
t342 = (t175 * t240 - t224 * t234) * t214;
t176 = mrSges(2,1) + t282;
t341 = (t176 * t242 - t224 * t236) * t214;
t223 = Ifges(3,1) - Ifges(3,2);
t308 = t231 * t237;
t162 = -t207 * t223 + t308 * t397 + t372;
t340 = t162 * t208;
t306 = t233 * t239;
t163 = -t209 * t223 + t306 * t397 + t372;
t339 = t163 * t210;
t304 = t235 * t241;
t164 = -t211 * t223 + t304 * t397 + t372;
t338 = t164 * t212;
t177 = Ifges(3,5) * t231 + Ifges(3,6) * t237;
t337 = t177 * t208;
t178 = Ifges(3,5) * t233 + Ifges(3,6) * t239;
t336 = t178 * t210;
t179 = Ifges(3,5) * t235 + Ifges(3,6) * t241;
t335 = t179 * t212;
t325 = t214 * t231;
t324 = t214 * t233;
t323 = t214 * t235;
t319 = t214 * t240;
t317 = t214 * t242;
t309 = t216 * t248;
t302 = pkin(2) * t358;
t299 = pkin(5) * t357;
t298 = pkin(5) * t355;
t297 = pkin(5) * t353;
t296 = t208 * t352;
t295 = t210 * t351;
t294 = t212 * t350;
t290 = t208 * t343;
t289 = t210 * t342;
t288 = t212 * t341;
t286 = t214 * t305;
t285 = t214 * t303;
t281 = Ifges(3,3) * t293;
t280 = Ifges(3,3) * t292;
t279 = Ifges(3,3) * t291;
t269 = mrSges(3,1) * t231 + mrSges(3,2) * t237;
t130 = -t269 * t214 * t232 + t284 * t216;
t278 = t130 * t293;
t268 = mrSges(3,1) * t233 + mrSges(3,2) * t239;
t131 = -t268 * t214 * t234 + t283 * t216;
t277 = t131 * t292;
t267 = mrSges(3,1) * t235 + mrSges(3,2) * t241;
t132 = -t267 * t214 * t236 + t282 * t216;
t276 = t132 * t291;
t275 = t177 * t293;
t274 = t178 * t292;
t273 = t179 * t291;
t272 = t293 * t367;
t271 = t292 * t366;
t270 = t291 * t365;
t266 = t133 * t214 + t216 * t88;
t265 = t134 * t214 + t216 * t89;
t264 = t135 * t214 + t216 * t90;
t263 = t103 * t213 + t104 * t215;
t262 = t105 * t213 + t106 * t215;
t261 = t107 * t213 + t108 * t215;
t257 = pkin(2) * t325 - t168 * t216;
t256 = pkin(2) * t324 - t169 * t216;
t255 = pkin(2) * t323 - t170 * t216;
t228 = xDDP(3);
t229 = xDDP(2);
t230 = xDDP(1);
t254 = t228 * t73 + t229 * t74 - t230 * t367;
t253 = t228 * t75 + t229 * t76 - t230 * t366;
t252 = t228 * t77 + t229 * t78 - t230 * t365;
t251 = t266 * t232 + t263 * t238;
t250 = t265 * t234 + t262 * t240;
t249 = t264 * t236 + t261 * t242;
t206 = m(1) + m(2) + m(3);
t173 = pkin(5) * t236 + t242 * t378;
t172 = pkin(5) * t234 + t240 * t379;
t171 = pkin(5) * t232 + t238 * t380;
t126 = -t173 * t213 + t215 * t255;
t125 = -t172 * t213 + t256 * t215;
t124 = -t171 * t213 + t257 * t215;
t123 = t173 * t215 + t213 * t255;
t122 = t172 * t215 + t256 * t213;
t121 = t171 * t215 + t257 * t213;
t96 = -t123 * t187 + t126 * t193;
t95 = -t122 * t186 + t125 * t192;
t94 = -t121 * t185 + t124 * t191;
t93 = (t152 * t255 + t155 * t173) * t205 + t202 * t146;
t92 = (t256 * t151 + t154 * t172) * t204 + t201 * t145;
t91 = (t257 * t150 + t153 * t171) * t203 + t200 * t144;
t87 = -t146 * t205 + (t123 * t193 + t126 * t187) * t202;
t86 = -t145 * t204 + (t122 * t192 + t125 * t186) * t201;
t85 = -t144 * t203 + (t121 * t191 + t124 * t185) * t200;
t72 = -t190 * t96 - t196 * t87;
t71 = t190 * t87 - t196 * t96;
t70 = -t189 * t95 - t195 * t86;
t69 = t189 * t86 - t195 * t95;
t68 = -t188 * t94 - t194 * t85;
t67 = t188 * t85 - t194 * t94;
t63 = -t179 * t270 + (-t164 * t294 + t93 * t341) * t394;
t62 = -t178 * t271 + (-t163 * t295 + t92 * t342) * t395;
t61 = -t177 * t272 + (-t162 * t296 + t91 * t343) * t141;
t60 = -t132 * t270 + (t206 * t93 - t288 * t350) * t394;
t59 = -t131 * t271 + (t206 * t92 - t289 * t351) * t395;
t58 = -t130 * t272 + (t206 * t91 - t290 * t352) * t141;
t54 = t57 ^ 2;
t53 = t56 ^ 2;
t52 = t55 ^ 2;
t48 = t51 ^ 2;
t47 = t50 ^ 2;
t46 = t49 ^ 2;
t39 = t77 * t273 + (t81 * t338 + t72 * t341) * t394;
t38 = t78 * t273 + (t84 * t338 + t71 * t341) * t394;
t37 = t75 * t274 + (t80 * t339 + t70 * t342) * t395;
t36 = t76 * t274 + (t83 * t339 + t69 * t342) * t395;
t35 = t73 * t275 + (t79 * t340 + t68 * t343) * t141;
t34 = t74 * t275 + (t82 * t340 + t67 * t343) * t141;
t33 = t77 * t276 + (t206 * t72 + t81 * t288) * t394;
t32 = t78 * t276 + (t206 * t71 + t84 * t288) * t394;
t31 = t75 * t277 + (t206 * t70 + t80 * t289) * t395;
t30 = t76 * t277 + (t206 * t69 + t83 * t289) * t395;
t29 = t73 * t278 + (t206 * t68 + t79 * t290) * t141;
t28 = t74 * t278 + (t206 * t67 + t82 * t290) * t141;
t24 = t297 - t386;
t23 = t298 - t387;
t22 = t299 - t388;
t19 = -pkin(5) * t302 + (t207 * t247 + t246) * t55;
t18 = (t24 * t386 - t373) * t394;
t17 = (t23 * t387 - t374) * t395;
t16 = (t19 * t55 - t22 * t388) * t396;
t15 = (t309 * t373 + (-t51 * t170 * t323 + t216 * (t51 * t382 - t297)) * t51) * t344;
t14 = (t309 * t374 + (-t50 * t169 * t324 + t216 * (t50 * t383 - t298)) * t50) * t345;
t13 = (t19 * t309 * t368 + (-t49 * t168 * t325 + t216 * (t49 * t384 - t299)) * t141 * t49) * t208;
t12 = (((t216 * t51 + t57 * t317) * t382 - (-pkin(5) * t57 + t300) * t285 + t216 * t24) * t57 - (-t51 * t317 + (-t211 * t216 + t235 * t285 + t216) * t57) * t386) * t344;
t11 = (((t216 * t50 + t56 * t319) * t383 - (-pkin(5) * t56 + t301) * t286 + t216 * t23) * t56 + (t50 * t319 + (t209 * t216 - t233 * t286 - t216) * t56) * t387) * t345;
t10 = (((t216 * t49 + t55 * t321) * t384 - (-pkin(5) * t55 + t302) * t287 + t216 * t22) * t368 + (t49 * t321 + (t207 * t216 - t231 * t287 - t216) * t55) * t396 * t388) * t208;
t9 = -t132 * t18 - t179 * t12 - Ifges(3,3) * t15 - t54 * (Ifges(3,4) * t391 + t223 * t304 - Ifges(3,4)) + (t400 * mrSges(3,1) + t249 * mrSges(3,2)) * t241 + (t249 * mrSges(3,1) - t400 * mrSges(3,2)) * t235;
t8 = -t131 * t17 - t178 * t11 - Ifges(3,3) * t14 - t53 * (Ifges(3,4) * t392 + t223 * t306 - Ifges(3,4)) + (t399 * mrSges(3,1) + t250 * mrSges(3,2)) * t239 + (t250 * mrSges(3,1) - t399 * mrSges(3,2)) * t233;
t7 = t130 * t16 - t177 * t10 - Ifges(3,3) * t13 - t52 * (Ifges(3,4) * t393 + t223 * t308 - Ifges(3,4)) + (t398 * mrSges(3,1) + t251 * mrSges(3,2)) * t237 + (t251 * mrSges(3,1) - t398 * mrSges(3,2)) * t231;
t6 = -t18 * t341 - t164 * t12 - t179 * t15 + 0.2e1 * t51 * ((t223 * t353 + t51 * t390) * t241 + t354 * t389 + (t391 - 0.1e1) * t57 * Ifges(3,4)) + (-t264 * t176 + t261 * t224) * t242 + (t261 * t176 + t264 * t224) * t236;
t5 = -t17 * t342 - t163 * t11 - t178 * t14 + 0.2e1 * t50 * ((t223 * t355 + t50 * t390) * t239 + t356 * t389 + (t392 - 0.1e1) * t56 * Ifges(3,4)) + (-t265 * t175 + t262 * t224) * t240 + (t262 * t175 + t265 * t224) * t234;
t4 = t16 * t343 - t162 * t10 - t177 * t13 + 0.2e1 * t49 * ((t223 * t357 + t49 * t390) * t237 + t358 * t389 + (t393 - 0.1e1) * t55 * Ifges(3,4)) + (-t266 * t174 + t263 * t224) * t238 + (t263 * t174 + t266 * t224) * t232;
t3 = -t12 * t341 - t132 * t15 + ((-mrSges(2,1) * t54 - t282 * (t54 + t48)) * t236 - 0.2e1 * t57 * t242 * (t267 * t51 + t57 * t385)) * t214 - t216 * t48 * t267 + (-t18 - t135) * t206;
t2 = -t11 * t342 - t131 * t14 + ((-mrSges(2,1) * t53 - t283 * (t53 + t47)) * t234 - 0.2e1 * t56 * t240 * (t268 * t50 + t56 * t385)) * t214 - t216 * t47 * t268 + (-t17 - t134) * t206;
t1 = -t10 * t343 - t130 * t13 + ((-mrSges(2,1) * t52 - t284 * (t52 + t46)) * t232 - 0.2e1 * t55 * t238 * (t269 * t49 + t55 * t385)) * t214 - t216 * t46 * t269 + (t16 - t133) * t206;
t20 = [(-g(1) + t230) * m(4) + ((-t63 * t294 + t60 * t93) * t230 + (t63 * t359 + t60 * t71) * t229 + (t63 * t360 + t60 * t72) * t228 + t93 * t3 - t6 * t294) * t394 + ((-t62 * t295 + t59 * t92) * t230 + (t62 * t361 + t59 * t69) * t229 + (t62 * t362 + t59 * t70) * t228 + t92 * t2 - t5 * t295) * t395 + ((-t9 * t365 + t252 * (-Ifges(3,3) * t270 + (t132 * t93 - t179 * t294) * t394)) * t344 + (-t8 * t366 + t253 * (-Ifges(3,3) * t271 + (t131 * t92 - t178 * t295) * t395)) * t345 + (-t7 * t367 + t254 * (-Ifges(3,3) * t272 + (t130 * t91 - t177 * t296) * t141)) * t346) * t248 + ((-t61 * t296 + t58 * t91) * t230 + (t61 * t363 + t58 * t67) * t229 + (t61 * t364 + t58 * t68) * t228 + t91 * t1 - t4 * t296) * t141; (-g(2) + t229) * m(4) + ((-t38 * t294 + t32 * t93) * t230 + (t32 * t71 + t38 * t359) * t229 + (t32 * t72 + t38 * t360) * t228 + t71 * t3 + t6 * t359) * t394 + ((-t36 * t295 + t30 * t92) * t230 + (t30 * t69 + t36 * t361) * t229 + (t30 * t70 + t36 * t362) * t228 + t69 * t2 + t5 * t361) * t395 + ((t78 * t9 + t252 * (t78 * t279 + (t132 * t71 + t84 * t335) * t394)) * t344 + (t76 * t8 + t253 * (t76 * t280 + (t131 * t69 + t83 * t336) * t395)) * t345 + (t7 * t74 + t254 * (t74 * t281 + (t130 * t67 + t82 * t337) * t141)) * t346) * t248 + ((t28 * t91 - t34 * t296) * t230 + (t28 * t67 + t34 * t363) * t229 + (t28 * t68 + t34 * t364) * t228 + t67 * t1 + t4 * t363) * t141; (-g(3) + t228) * m(4) + ((-t39 * t294 + t33 * t93) * t230 + (t33 * t71 + t39 * t359) * t229 + (t33 * t72 + t39 * t360) * t228 + t72 * t3 + t6 * t360) * t394 + ((-t37 * t295 + t31 * t92) * t230 + (t31 * t69 + t37 * t361) * t229 + (t31 * t70 + t37 * t362) * t228 + t70 * t2 + t5 * t362) * t395 + ((t77 * t9 + t252 * (t77 * t279 + (t132 * t72 + t81 * t335) * t394)) * t344 + (t75 * t8 + t253 * (t75 * t280 + (t131 * t70 + t80 * t336) * t395)) * t345 + (t7 * t73 + t254 * (t73 * t281 + (t130 * t68 + t79 * t337) * t141)) * t346) * t248 + ((t29 * t91 - t35 * t296) * t230 + (t29 * t67 + t35 * t363) * t229 + (t29 * t68 + t35 * t364) * t228 + t68 * t1 + t4 * t364) * t141;];
tauX  = t20;
