% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:25
% EndTime: 2020-08-06 18:21:26
% DurationCPUTime: 0.36s
% Computational Cost: add. (1203->136), mult. (1239->245), div. (189->7), fcn. (945->18), ass. (0->97)
t347 = 2 * mrSges(3,1);
t346 = 2 * mrSges(3,2);
t345 = 2 * mrSges(2,3);
t344 = -2 * Ifges(3,4);
t331 = (pkin(1) + pkin(5));
t343 = (mrSges(3,3) - mrSges(2,2));
t334 = 0.1e1 / pkin(3);
t342 = Ifges(3,3) * t334;
t319 = sin(qJ(3,3));
t298 = t319 * pkin(3) + qJ(2,3);
t295 = 0.1e1 / t298;
t312 = 0.1e1 / t319;
t341 = t295 * t312;
t321 = sin(qJ(3,2));
t299 = t321 * pkin(3) + qJ(2,2);
t296 = 0.1e1 / t299;
t313 = 0.1e1 / t321;
t340 = t296 * t313;
t323 = sin(qJ(3,1));
t300 = t323 * pkin(3) + qJ(2,1);
t297 = 0.1e1 / t300;
t314 = 0.1e1 / t323;
t339 = t297 * t314;
t325 = cos(qJ(3,3));
t338 = t312 * t325;
t327 = cos(qJ(3,2));
t337 = t313 * t327;
t329 = cos(qJ(3,1));
t336 = t314 * t329;
t332 = (m(2) + m(3));
t335 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 2 * (m(3) * pkin(5) + t343) * pkin(1) + t332 * pkin(1) ^ 2 + 2 * mrSges(3,3) * pkin(5) + m(3) * pkin(5) ^ 2;
t330 = cos(qJ(1,1));
t328 = cos(qJ(1,2));
t326 = cos(qJ(1,3));
t324 = sin(qJ(1,1));
t322 = sin(qJ(1,2));
t320 = sin(qJ(1,3));
t318 = Ifges(3,1) - Ifges(3,2);
t317 = legFrame(1,3);
t316 = legFrame(2,3);
t315 = legFrame(3,3);
t310 = pkin(6) + t331;
t308 = cos(t317);
t307 = cos(t316);
t306 = cos(t315);
t305 = sin(t317);
t304 = sin(t316);
t303 = sin(t315);
t302 = -t331 * mrSges(3,1) + Ifges(3,5);
t301 = t331 * mrSges(3,2) - Ifges(3,6);
t293 = -m(2) * pkin(1) - t331 * m(3) - t343;
t292 = t329 * mrSges(3,1) - t323 * mrSges(3,2);
t291 = t327 * mrSges(3,1) - t321 * mrSges(3,2);
t290 = t325 * mrSges(3,1) - t319 * mrSges(3,2);
t289 = t305 * t330 + t308 * t324;
t288 = t304 * t328 + t307 * t322;
t287 = t303 * t326 + t306 * t320;
t286 = -t305 * t324 + t308 * t330;
t285 = -t304 * t322 + t307 * t328;
t284 = -t303 * t320 + t306 * t326;
t283 = -t300 * t330 + t310 * t324;
t282 = -t299 * t328 + t310 * t322;
t281 = -t298 * t326 + t310 * t320;
t280 = t324 * t300 + t310 * t330;
t279 = t322 * t299 + t310 * t328;
t278 = t320 * t298 + t310 * t326;
t277 = t323 * t301 + t302 * t329;
t276 = t321 * t301 + t302 * t327;
t275 = t319 * t301 + t302 * t325;
t274 = (-t292 * t334 + t329 * t332) * t314;
t273 = (-t291 * t334 + t327 * t332) * t313;
t272 = (-t290 * t334 + t325 * t332) * t312;
t271 = (-t277 * t334 + t293 * t329) * t314;
t270 = (-t276 * t334 + t293 * t327) * t313;
t269 = (-t275 * t334 + t293 * t325) * t312;
t268 = t305 * t280 + t283 * t308;
t267 = t304 * t279 + t282 * t307;
t266 = t303 * t278 + t281 * t306;
t265 = t280 * t308 - t305 * t283;
t264 = t279 * t307 - t304 * t282;
t263 = t278 * t306 - t303 * t281;
t262 = (t332 * qJ(2,1) + t323 * t347 + t345) * qJ(2,1) + (qJ(2,1) * t346 + t318 * t329 + t323 * t344) * t329 + t335;
t261 = (t332 * qJ(2,2) + t321 * t347 + t345) * qJ(2,2) + (qJ(2,2) * t346 + t318 * t327 + t321 * t344) * t327 + t335;
t260 = (t332 * qJ(2,3) + t319 * t347 + t345) * qJ(2,3) + (qJ(2,3) * t346 + t318 * t325 + t319 * t344) * t325 + t335;
t259 = (t268 * t332 + t289 * t293) * t297;
t258 = (t267 * t332 + t288 * t293) * t296;
t257 = (t266 * t332 + t287 * t293) * t295;
t256 = (t265 * t332 + t286 * t293) * t297;
t255 = (t264 * t332 + t285 * t293) * t296;
t254 = (t263 * t332 + t284 * t293) * t295;
t253 = (t262 * t289 + t268 * t293) * t297;
t252 = (t261 * t288 + t267 * t293) * t296;
t251 = (t260 * t287 + t266 * t293) * t295;
t250 = (t262 * t286 + t265 * t293) * t297;
t249 = (t261 * t285 + t264 * t293) * t296;
t248 = (t260 * t284 + t263 * t293) * t295;
t1 = [m(4) + (t250 * t286 + t256 * t265) * t297 + (t249 * t285 + t255 * t264) * t296 + (t248 * t284 + t254 * t263) * t295, (t250 * t289 + t256 * t268) * t297 + (t249 * t288 + t255 * t267) * t296 + (t248 * t287 + t254 * t266) * t295, t254 * t338 + t255 * t337 + t256 * t336 + (-(t265 * t292 + t277 * t286) * t339 - (t264 * t291 + t276 * t285) * t340 - (t263 * t290 + t275 * t284) * t341) * t334; (t253 * t286 + t259 * t265) * t297 + (t252 * t285 + t258 * t264) * t296 + (t251 * t284 + t257 * t263) * t295, m(4) + (t253 * t289 + t259 * t268) * t297 + (t252 * t288 + t258 * t267) * t296 + (t251 * t287 + t257 * t266) * t295, t257 * t338 + t258 * t337 + t259 * t336 + (-(t268 * t292 + t277 * t289) * t339 - (t267 * t291 + t276 * t288) * t340 - (t266 * t290 + t275 * t287) * t341) * t334; (t265 * t274 + t271 * t286) * t297 + (t264 * t273 + t270 * t285) * t296 + (t263 * t272 + t269 * t284) * t295, (t268 * t274 + t271 * t289) * t297 + (t267 * t273 + t270 * t288) * t296 + (t266 * t272 + t269 * t287) * t295, m(4) + (t274 * t329 - (t292 * t329 - t342) * t334 * t314) * t314 + (t273 * t327 - (t291 * t327 - t342) * t334 * t313) * t313 + (t272 * t325 - (t290 * t325 - t342) * t334 * t312) * t312;];
MX  = t1;
