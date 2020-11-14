% Calculate Gravitation load for parallel robot
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRR1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:37
% EndTime: 2019-05-03 14:47:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (132->55), mult. (239->104), div. (36->4), fcn. (218->14), ass. (0->56)
t357 = legFrame(3,3);
t346 = sin(t357);
t349 = cos(t357);
t340 = -t346 * g(1) + t349 * g(2);
t343 = t349 * g(1) + t346 * g(2);
t360 = sin(qJ(2,3));
t354 = 0.1e1 / t360;
t363 = cos(qJ(2,3));
t382 = ((mrSges(2,1) * t343 + mrSges(2,2) * t340) * t360 - t363 * (mrSges(2,1) * t340 - mrSges(2,2) * t343)) * t354;
t358 = legFrame(2,3);
t347 = sin(t358);
t350 = cos(t358);
t341 = -t347 * g(1) + t350 * g(2);
t344 = t350 * g(1) + t347 * g(2);
t361 = sin(qJ(2,2));
t355 = 0.1e1 / t361;
t364 = cos(qJ(2,2));
t381 = ((mrSges(2,1) * t344 + mrSges(2,2) * t341) * t361 - t364 * (mrSges(2,1) * t341 - mrSges(2,2) * t344)) * t355;
t359 = legFrame(1,3);
t348 = sin(t359);
t351 = cos(t359);
t342 = -t348 * g(1) + t351 * g(2);
t345 = t351 * g(1) + t348 * g(2);
t362 = sin(qJ(2,1));
t356 = 0.1e1 / t362;
t365 = cos(qJ(2,1));
t380 = ((mrSges(2,1) * t345 + mrSges(2,2) * t342) * t362 - t365 * (mrSges(2,1) * t342 - mrSges(2,2) * t345)) * t356;
t379 = t340 * t354;
t378 = t341 * t355;
t377 = t342 * t356;
t376 = 0.1e1 / pkin(2);
t375 = koppelP(1,1);
t374 = koppelP(2,1);
t373 = koppelP(3,1);
t372 = koppelP(1,2);
t371 = koppelP(2,2);
t370 = koppelP(3,2);
t369 = mrSges(3,1);
t368 = mrSges(3,2);
t367 = xP(3);
t366 = m(1) + m(2);
t353 = cos(t367);
t352 = sin(t367);
t339 = -t352 * t372 + t353 * t375;
t338 = -t352 * t371 + t353 * t374;
t337 = -t352 * t370 + t353 * t373;
t336 = -t352 * t375 - t353 * t372;
t335 = -t352 * t374 - t353 * t371;
t334 = -t352 * t373 - t353 * t370;
t333 = t348 * t365 + t362 * t351;
t332 = -t348 * t362 + t351 * t365;
t331 = t347 * t364 + t361 * t350;
t330 = -t347 * t361 + t350 * t364;
t329 = t346 * t363 + t360 * t349;
t328 = -t346 * t360 + t349 * t363;
t1 = [-g(1) * m(3) + (-t349 * t382 - t350 * t381 - t351 * t380) * t376 + (-t328 * t379 - t330 * t378 - t332 * t377) * t366; -g(2) * m(3) + (-t346 * t382 - t347 * t381 - t348 * t380) * t376 + (-t329 * t379 - t331 * t378 - t333 * t377) * t366; -(-g(1) * t369 - g(2) * t368) * t352 + t353 * (g(1) * t368 - g(2) * t369) + (-(t332 * t336 + t333 * t339) * t377 - (t330 * t335 + t331 * t338) * t378 - (t328 * t334 + t329 * t337) * t379) * t366 + ((-t336 * t351 - t339 * t348) * t380 + (-t335 * t350 - t338 * t347) * t381 + (-t334 * t349 - t337 * t346) * t382) * t376;];
taugX  = t1;
