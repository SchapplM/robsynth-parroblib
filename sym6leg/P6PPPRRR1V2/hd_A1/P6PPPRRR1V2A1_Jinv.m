% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d4,theta1,theta2,theta3]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-16 19:46
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PPPRRR1V2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPPRRR1V2A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPPRRR1V2A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPPRRR1V2A1_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPPRRR1V2A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPPRRR1V2A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 19:45:40
% EndTime: 2019-05-16 19:45:43
% DurationCPUTime: 3.90s
% Computational Cost: add. (1122->372), mult. (2976->800), div. (72->1), fcn. (3282->48), ass. (0->331)
t227 = xP(4);
t195 = sin(t227);
t226 = xP(5);
t197 = cos(t226);
t228 = koppelP(6,3);
t314 = t197 * t228;
t225 = xP(6);
t193 = sin(t225);
t194 = sin(t226);
t196 = cos(t225);
t198 = cos(t227);
t240 = koppelP(6,1);
t302 = t198 * t240;
t234 = koppelP(6,2);
t308 = t198 * t234;
t320 = t195 * t240;
t326 = t195 * t234;
t49 = (t194 * t320 + t308) * t196 + t193 * (-t194 * t326 + t302);
t374 = t195 * t314 - t49;
t229 = koppelP(5,3);
t313 = t197 * t229;
t241 = koppelP(5,1);
t301 = t198 * t241;
t235 = koppelP(5,2);
t307 = t198 * t235;
t319 = t195 * t241;
t325 = t195 * t235;
t50 = (t194 * t319 + t307) * t196 + t193 * (-t194 * t325 + t301);
t373 = t195 * t313 - t50;
t230 = koppelP(4,3);
t312 = t197 * t230;
t242 = koppelP(4,1);
t300 = t198 * t242;
t236 = koppelP(4,2);
t306 = t198 * t236;
t318 = t195 * t242;
t324 = t195 * t236;
t51 = (t194 * t318 + t306) * t196 + t193 * (-t194 * t324 + t300);
t372 = t195 * t312 - t51;
t231 = koppelP(3,3);
t311 = t197 * t231;
t243 = koppelP(3,1);
t299 = t198 * t243;
t237 = koppelP(3,2);
t305 = t198 * t237;
t317 = t195 * t243;
t323 = t195 * t237;
t52 = (t194 * t317 + t305) * t196 + t193 * (-t194 * t323 + t299);
t371 = t195 * t311 - t52;
t232 = koppelP(2,3);
t310 = t197 * t232;
t244 = koppelP(2,1);
t298 = t198 * t244;
t238 = koppelP(2,2);
t304 = t198 * t238;
t316 = t195 * t244;
t322 = t195 * t238;
t53 = (t194 * t316 + t304) * t196 + t193 * (-t194 * t322 + t298);
t370 = t195 * t310 - t53;
t233 = koppelP(1,3);
t309 = t197 * t233;
t245 = koppelP(1,1);
t297 = t198 * t245;
t239 = koppelP(1,2);
t303 = t198 * t239;
t315 = t195 * t245;
t321 = t195 * t239;
t54 = (t194 * t315 + t303) * t196 + t193 * (-t194 * t321 + t297);
t369 = t195 * t309 - t54;
t207 = legFrame(6,3);
t157 = sin(t207);
t169 = cos(t207);
t202 = sin(pkin(8));
t205 = cos(pkin(8));
t133 = t205 * t157 + t202 * t169;
t219 = legFrame(6,2);
t187 = cos(t219);
t262 = t187 * t133;
t208 = legFrame(5,3);
t158 = sin(t208);
t170 = cos(t208);
t134 = t205 * t158 + t202 * t170;
t220 = legFrame(5,2);
t188 = cos(t220);
t260 = t188 * t134;
t209 = legFrame(4,3);
t159 = sin(t209);
t171 = cos(t209);
t135 = t205 * t159 + t202 * t171;
t221 = legFrame(4,2);
t189 = cos(t221);
t258 = t189 * t135;
t210 = legFrame(3,3);
t160 = sin(t210);
t172 = cos(t210);
t136 = t205 * t160 + t202 * t172;
t222 = legFrame(3,2);
t190 = cos(t222);
t256 = t190 * t136;
t211 = legFrame(2,3);
t161 = sin(t211);
t173 = cos(t211);
t137 = t205 * t161 + t202 * t173;
t223 = legFrame(2,2);
t191 = cos(t223);
t254 = t191 * t137;
t212 = legFrame(1,3);
t162 = sin(t212);
t174 = cos(t212);
t138 = t205 * t162 + t202 * t174;
t224 = legFrame(1,2);
t192 = cos(t224);
t252 = t192 * t138;
t201 = sin(pkin(9));
t203 = sin(pkin(4));
t295 = t201 * t203;
t127 = t157 * t202 - t205 * t169;
t181 = sin(t219);
t368 = t127 * t181;
t367 = t127 * t187;
t206 = cos(pkin(4));
t366 = t127 * t206;
t128 = t158 * t202 - t205 * t170;
t182 = sin(t220);
t365 = t128 * t182;
t364 = t128 * t188;
t363 = t128 * t206;
t129 = t159 * t202 - t205 * t171;
t183 = sin(t221);
t362 = t129 * t183;
t361 = t129 * t189;
t360 = t129 * t206;
t130 = t160 * t202 - t205 * t172;
t184 = sin(t222);
t359 = t130 * t184;
t358 = t130 * t190;
t357 = t130 * t206;
t131 = t161 * t202 - t205 * t173;
t185 = sin(t223);
t356 = t131 * t185;
t355 = t131 * t191;
t354 = t131 * t206;
t132 = t162 * t202 - t205 * t174;
t186 = sin(t224);
t353 = t132 * t186;
t352 = t132 * t192;
t351 = t132 * t206;
t269 = t193 * t234 - t196 * t240;
t350 = t269 * t197;
t268 = t193 * t235 - t196 * t241;
t349 = t268 * t197;
t267 = t193 * t236 - t196 * t242;
t348 = t267 * t197;
t266 = t193 * t237 - t196 * t243;
t347 = t266 * t197;
t265 = t193 * t238 - t196 * t244;
t346 = t265 * t197;
t264 = t193 * t239 - t196 * t245;
t345 = t264 * t197;
t344 = t181 * t202;
t343 = t181 * t205;
t342 = t182 * t202;
t341 = t182 * t205;
t340 = t183 * t202;
t339 = t183 * t205;
t338 = t184 * t202;
t337 = t184 * t205;
t336 = t185 * t202;
t335 = t185 * t205;
t334 = t186 * t202;
t333 = t186 * t205;
t332 = t187 * t203;
t331 = t188 * t203;
t330 = t189 * t203;
t329 = t190 * t203;
t328 = t191 * t203;
t327 = t192 * t203;
t296 = 0.1e1 / t295;
t294 = t201 * t206;
t55 = -(t194 * t308 + t320) * t193 + (t194 * t302 - t326) * t196;
t287 = t198 * t314 - t55;
t56 = -(t194 * t307 + t319) * t193 + (t194 * t301 - t325) * t196;
t286 = t198 * t313 - t56;
t57 = -(t194 * t306 + t318) * t193 + (t194 * t300 - t324) * t196;
t285 = t198 * t312 - t57;
t58 = -(t194 * t305 + t317) * t193 + (t194 * t299 - t323) * t196;
t284 = t198 * t311 - t58;
t59 = -(t194 * t304 + t316) * t193 + (t194 * t298 - t322) * t196;
t283 = t198 * t310 - t59;
t60 = -(t194 * t303 + t315) * t193 + (t194 * t297 - t321) * t196;
t282 = t198 * t309 - t60;
t213 = legFrame(6,1);
t163 = sin(t213);
t175 = cos(t213);
t103 = t163 * t343 + t175 * t202;
t91 = -t163 * t344 + t175 * t205;
t281 = -t103 * t157 + t169 * t91;
t104 = -t163 * t202 + t175 * t343;
t92 = t163 * t205 + t175 * t344;
t280 = -t104 * t157 - t169 * t92;
t214 = legFrame(5,1);
t164 = sin(t214);
t176 = cos(t214);
t105 = t164 * t341 + t176 * t202;
t93 = -t164 * t342 + t176 * t205;
t279 = -t105 * t158 + t170 * t93;
t106 = -t164 * t202 + t176 * t341;
t94 = t164 * t205 + t176 * t342;
t278 = -t106 * t158 - t170 * t94;
t215 = legFrame(4,1);
t165 = sin(t215);
t177 = cos(t215);
t107 = t165 * t339 + t177 * t202;
t95 = -t165 * t340 + t177 * t205;
t277 = -t107 * t159 + t171 * t95;
t108 = -t165 * t202 + t177 * t339;
t96 = t165 * t205 + t177 * t340;
t276 = -t108 * t159 - t171 * t96;
t216 = legFrame(3,1);
t166 = sin(t216);
t178 = cos(t216);
t109 = t166 * t337 + t178 * t202;
t97 = -t166 * t338 + t178 * t205;
t275 = -t109 * t160 + t172 * t97;
t110 = -t166 * t202 + t178 * t337;
t98 = t166 * t205 + t178 * t338;
t274 = -t110 * t160 - t172 * t98;
t217 = legFrame(2,1);
t167 = sin(t217);
t179 = cos(t217);
t111 = t167 * t335 + t179 * t202;
t99 = -t167 * t336 + t179 * t205;
t273 = -t111 * t161 + t173 * t99;
t100 = t167 * t205 + t179 * t336;
t112 = -t167 * t202 + t179 * t335;
t272 = -t100 * t173 - t112 * t161;
t218 = legFrame(1,1);
t168 = sin(t218);
t180 = cos(t218);
t101 = -t168 * t334 + t180 * t205;
t113 = t168 * t333 + t180 * t202;
t271 = t101 * t174 - t113 * t162;
t102 = t168 * t205 + t180 * t334;
t114 = -t168 * t202 + t180 * t333;
t270 = -t102 * t174 - t114 * t162;
t73 = t194 * t228 - t350;
t74 = t194 * t229 - t349;
t75 = t194 * t230 - t348;
t76 = t194 * t231 - t347;
t77 = t194 * t232 - t346;
t78 = t194 * t233 - t345;
t263 = t187 * t73;
t261 = t188 * t74;
t259 = t189 * t75;
t257 = t190 * t76;
t255 = t191 * t77;
t253 = t192 * t78;
t251 = t197 * t262;
t250 = t197 * t260;
t249 = t197 * t258;
t248 = t197 * t256;
t247 = t197 * t254;
t246 = t197 * t252;
t204 = cos(pkin(9));
t150 = t193 * t245 + t196 * t239;
t149 = t193 * t244 + t196 * t238;
t148 = t193 * t243 + t196 * t237;
t147 = t193 * t242 + t196 * t236;
t146 = t193 * t241 + t196 * t235;
t145 = t193 * t240 + t196 * t234;
t72 = t264 * t194 + t309;
t71 = t265 * t194 + t310;
t70 = t266 * t194 + t311;
t69 = t267 * t194 + t312;
t68 = t268 * t194 + t313;
t67 = t269 * t194 + t314;
t66 = t138 * t206 * t186 + t327;
t65 = t137 * t206 * t185 + t328;
t64 = t136 * t206 * t184 + t329;
t63 = t135 * t206 * t183 + t330;
t62 = t134 * t206 * t182 + t331;
t61 = t133 * t206 * t181 + t332;
t48 = t150 * t333 + t202 * t72;
t47 = t149 * t335 + t202 * t71;
t46 = t148 * t337 + t202 * t70;
t45 = t147 * t339 + t202 * t69;
t44 = t146 * t341 + t202 * t68;
t43 = t145 * t343 + t202 * t67;
t42 = -t150 * t202 + t72 * t333;
t41 = -t149 * t202 + t71 * t335;
t40 = -t148 * t202 + t70 * t337;
t39 = -t147 * t202 + t69 * t339;
t38 = -t146 * t202 + t68 * t341;
t37 = -t145 * t202 + t67 * t343;
t36 = -t150 * t334 + t205 * t72;
t35 = -t149 * t336 + t205 * t71;
t34 = -t148 * t338 + t205 * t70;
t33 = -t147 * t340 + t205 * t69;
t32 = -t146 * t342 + t205 * t68;
t31 = -t145 * t344 + t205 * t67;
t30 = t150 * t205 + t72 * t334;
t29 = t149 * t205 + t71 * t336;
t28 = t148 * t205 + t70 * t338;
t27 = t147 * t205 + t69 * t340;
t26 = t146 * t205 + t68 * t342;
t25 = t145 * t205 + t67 * t344;
t24 = t162 * t36 + t174 * t48;
t23 = t161 * t35 + t173 * t47;
t22 = t160 * t34 + t172 * t46;
t21 = t159 * t33 + t171 * t45;
t20 = t158 * t32 + t170 * t44;
t19 = t157 * t31 + t169 * t43;
t18 = -t162 * t30 + t174 * t42;
t17 = -t161 * t29 + t173 * t41;
t16 = -t160 * t28 + t172 * t40;
t15 = -t159 * t27 + t171 * t39;
t14 = -t158 * t26 + t170 * t38;
t13 = -t157 * t25 + t169 * t37;
t12 = -t150 * t327 + (-t162 * t48 + t174 * t36) * t206;
t11 = -t149 * t328 + (-t161 * t47 + t173 * t35) * t206;
t10 = -t148 * t329 + (-t160 * t46 + t172 * t34) * t206;
t9 = -t147 * t330 + (-t159 * t45 + t171 * t33) * t206;
t8 = -t146 * t331 + (-t158 * t44 + t170 * t32) * t206;
t7 = -t145 * t332 + (-t157 * t43 + t169 * t31) * t206;
t6 = t72 * t327 + (t162 * t42 + t174 * t30) * t206;
t5 = t71 * t328 + (t161 * t41 + t173 * t29) * t206;
t4 = t70 * t329 + (t160 * t40 + t172 * t28) * t206;
t3 = t69 * t330 + (t159 * t39 + t171 * t27) * t206;
t2 = t68 * t331 + (t158 * t38 + t170 * t26) * t206;
t1 = t67 * t332 + (t157 * t37 + t169 * t25) * t206;
t79 = [((-t132 * t204 - t138 * t294) * t192 + t186 * t295) * t296, ((-t168 * t66 - t180 * t351) * t201 - t204 * (-t138 * t180 + t168 * t353)) * t296, ((-t168 * t351 + t180 * t66) * t201 + t204 * (t138 * t168 + t180 * t353)) * t296, (((-t12 * t198 - t195 * t6) * t180 + t168 * (-t12 * t195 + t198 * t6)) * t201 - ((-t18 * t195 + t198 * t24) * t180 + t168 * (t18 * t198 + t195 * t24)) * t204) * t296, ((-t282 * t352 + t78 * (-t102 * t162 + t114 * t174)) * t204 + (-t180 * t253 + t186 * t282) * t295 + (t60 * t252 - t270 * t345 + (t194 * t270 - t198 * t246) * t233) * t294) * t296, (-(t369 * t352 - t78 * (t101 * t162 + t113 * t174)) * t204 + (-t168 * t253 + t369 * t186) * t295 + (t54 * t252 - t271 * t345 + (t194 * t271 - t195 * t246) * t233) * t294) * t296; ((-t131 * t204 - t137 * t294) * t191 + t185 * t295) * t296, ((-t167 * t65 - t179 * t354) * t201 - t204 * (-t137 * t179 + t167 * t356)) * t296, ((-t167 * t354 + t179 * t65) * t201 + t204 * (t137 * t167 + t179 * t356)) * t296, (((-t11 * t198 - t195 * t5) * t179 + t167 * (-t11 * t195 + t198 * t5)) * t201 - ((-t17 * t195 + t198 * t23) * t179 + t167 * (t17 * t198 + t195 * t23)) * t204) * t296, ((-t283 * t355 + t77 * (-t100 * t161 + t112 * t173)) * t204 + (-t179 * t255 + t185 * t283) * t295 + (t59 * t254 - t272 * t346 + (t194 * t272 - t198 * t247) * t232) * t294) * t296, (-(t370 * t355 - t77 * (t111 * t173 + t161 * t99)) * t204 + (-t167 * t255 + t370 * t185) * t295 + (t53 * t254 - t273 * t346 + (t194 * t273 - t195 * t247) * t232) * t294) * t296; ((-t130 * t204 - t136 * t294) * t190 + t184 * t295) * t296, ((-t166 * t64 - t178 * t357) * t201 - t204 * (-t136 * t178 + t166 * t359)) * t296, ((-t166 * t357 + t178 * t64) * t201 + t204 * (t136 * t166 + t178 * t359)) * t296, (((-t10 * t198 - t195 * t4) * t178 + t166 * (-t10 * t195 + t198 * t4)) * t201 - ((-t16 * t195 + t198 * t22) * t178 + t166 * (t16 * t198 + t195 * t22)) * t204) * t296, ((-t284 * t358 + t76 * (t110 * t172 - t160 * t98)) * t204 + (-t178 * t257 + t184 * t284) * t295 + (t58 * t256 - t274 * t347 + (t194 * t274 - t198 * t248) * t231) * t294) * t296, (-(t371 * t358 - t76 * (t109 * t172 + t160 * t97)) * t204 + (-t166 * t257 + t371 * t184) * t295 + (t52 * t256 - t275 * t347 + (t194 * t275 - t195 * t248) * t231) * t294) * t296; ((-t129 * t204 - t135 * t294) * t189 + t183 * t295) * t296, ((-t165 * t63 - t177 * t360) * t201 - t204 * (-t135 * t177 + t165 * t362)) * t296, ((-t165 * t360 + t177 * t63) * t201 + t204 * (t135 * t165 + t177 * t362)) * t296, (((-t195 * t3 - t198 * t9) * t177 + t165 * (-t195 * t9 + t198 * t3)) * t201 - ((-t15 * t195 + t198 * t21) * t177 + t165 * (t15 * t198 + t195 * t21)) * t204) * t296, ((-t285 * t361 + t75 * (t108 * t171 - t159 * t96)) * t204 + (-t177 * t259 + t183 * t285) * t295 + (t57 * t258 - t276 * t348 + (t194 * t276 - t198 * t249) * t230) * t294) * t296, (-(t372 * t361 - t75 * (t107 * t171 + t159 * t95)) * t204 + (-t165 * t259 + t372 * t183) * t295 + (t51 * t258 - t277 * t348 + (t194 * t277 - t195 * t249) * t230) * t294) * t296; ((-t128 * t204 - t134 * t294) * t188 + t182 * t295) * t296, ((-t164 * t62 - t176 * t363) * t201 - t204 * (-t134 * t176 + t164 * t365)) * t296, ((-t164 * t363 + t176 * t62) * t201 + t204 * (t134 * t164 + t176 * t365)) * t296, (((-t195 * t2 - t198 * t8) * t176 + t164 * (-t195 * t8 + t198 * t2)) * t201 - ((-t14 * t195 + t198 * t20) * t176 + t164 * (t14 * t198 + t195 * t20)) * t204) * t296, ((-t286 * t364 + t74 * (t106 * t170 - t158 * t94)) * t204 + (-t176 * t261 + t286 * t182) * t295 + (t56 * t260 - t278 * t349 + (t194 * t278 - t198 * t250) * t229) * t294) * t296, (-(t373 * t364 - t74 * (t105 * t170 + t158 * t93)) * t204 + (-t164 * t261 + t373 * t182) * t295 + (t50 * t260 - t279 * t349 + (t194 * t279 - t195 * t250) * t229) * t294) * t296; ((-t127 * t204 - t133 * t294) * t187 + t181 * t295) * t296, ((-t163 * t61 - t175 * t366) * t201 - t204 * (-t133 * t175 + t163 * t368)) * t296, ((-t163 * t366 + t175 * t61) * t201 + t204 * (t133 * t163 + t175 * t368)) * t296, (((-t1 * t195 - t198 * t7) * t175 + t163 * (t1 * t198 - t195 * t7)) * t201 - ((-t13 * t195 + t19 * t198) * t175 + t163 * (t13 * t198 + t19 * t195)) * t204) * t296, ((-t287 * t367 + t73 * (t104 * t169 - t157 * t92)) * t204 + (-t175 * t263 + t181 * t287) * t295 + (t55 * t262 - t280 * t350 + (t194 * t280 - t198 * t251) * t228) * t294) * t296, (-(t374 * t367 - t73 * (t103 * t169 + t157 * t91)) * t204 + (-t163 * t263 + t374 * t181) * t295 + (t49 * t262 - t281 * t350 + (t194 * t281 - t195 * t251) * t228) * t294) * t296;];
Jinv  = t79;
