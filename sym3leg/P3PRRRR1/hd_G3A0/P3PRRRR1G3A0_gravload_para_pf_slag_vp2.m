% Calculate Gravitation load for parallel robot
% P3PRRRR1G3P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:22
% EndTime: 2020-03-09 21:02:22
% DurationCPUTime: 0.19s
% Computational Cost: add. (234->62), mult. (396->116), div. (45->10), fcn. (330->18), ass. (0->57)
t528 = mrSges(3,1) * g(3);
t496 = legFrame(3,2);
t482 = sin(t496);
t485 = cos(t496);
t476 = t482 * g(1) + t485 * g(2);
t500 = sin(qJ(2,3));
t489 = 0.1e1 / t500;
t527 = t476 * t489;
t497 = legFrame(2,2);
t483 = sin(t497);
t486 = cos(t497);
t477 = t483 * g(1) + t486 * g(2);
t502 = sin(qJ(2,2));
t490 = 0.1e1 / t502;
t526 = t477 * t490;
t498 = legFrame(1,2);
t484 = sin(t498);
t487 = cos(t498);
t478 = t484 * g(1) + t487 * g(2);
t504 = sin(qJ(2,1));
t491 = 0.1e1 / t504;
t525 = t478 * t491;
t499 = sin(qJ(3,3));
t524 = t489 * t499;
t501 = sin(qJ(3,2));
t523 = t490 * t501;
t503 = sin(qJ(3,1));
t522 = t491 * t503;
t479 = t485 * g(1) - t482 * g(2);
t495 = mrSges(2,2) - mrSges(3,3);
t506 = cos(qJ(2,3));
t505 = cos(qJ(3,3));
t515 = mrSges(3,1) * t505 - mrSges(3,2) * t499 + mrSges(2,1);
t473 = (t495 * t506 + t515 * t500) * t479 + (t495 * t500 - t515 * t506) * t476;
t492 = 0.1e1 / t505;
t521 = t473 * t489 * t492;
t480 = t486 * g(1) - t483 * g(2);
t508 = cos(qJ(2,2));
t507 = cos(qJ(3,2));
t514 = mrSges(3,1) * t507 - mrSges(3,2) * t501 + mrSges(2,1);
t474 = (t495 * t508 + t514 * t502) * t480 + (t495 * t502 - t514 * t508) * t477;
t493 = 0.1e1 / t507;
t520 = t474 * t490 * t493;
t481 = t487 * g(1) - t484 * g(2);
t510 = cos(qJ(2,1));
t509 = cos(qJ(3,1));
t513 = mrSges(3,1) * t509 - mrSges(3,2) * t503 + mrSges(2,1);
t475 = (t495 * t510 + t513 * t504) * t481 + (t495 * t504 - t513 * t510) * t478;
t494 = 0.1e1 / t509;
t519 = t475 * t491 * t494;
t518 = t476 * t500 + t479 * t506;
t517 = t477 * t502 + t480 * t508;
t516 = t478 * t504 + t481 * t510;
t512 = 0.1e1 / pkin(2);
t511 = mrSges(3,2) * g(3);
t488 = m(1) + m(2) + m(3);
t1 = [-g(1) * m(4) + (-t485 * t521 - t486 * t520 - t487 * t519) * t512 + (-(t504 * t484 + t487 * t510) * t525 - (t502 * t483 + t486 * t508) * t526 - (t500 * t482 + t485 * t506) * t527) * t488; -g(2) * m(4) + (t482 * t521 + t483 * t520 + t484 * t519) * t512 + (-(-t484 * t510 + t487 * t504) * t525 - (-t483 * t508 + t486 * t502) * t526 - (-t482 * t506 + t485 * t500) * t527) * t488; -g(3) * m(4) + (-t492 * t476 * t524 - t493 * t477 * t523 - t494 * t478 * t522) * t488 + (-t510 / t509 ^ 2 * t475 * t522 + t494 * ((t516 * mrSges(3,2) - t528) * t509 + t503 * (t516 * mrSges(3,1) + t511)) - t508 / t507 ^ 2 * t474 * t523 + t493 * ((t517 * mrSges(3,2) - t528) * t507 + t501 * (t517 * mrSges(3,1) + t511)) - t506 / t505 ^ 2 * t473 * t524 + t492 * ((t518 * mrSges(3,2) - t528) * t505 + t499 * (t518 * mrSges(3,1) + t511))) * t512;];
taugX  = t1;
