% Calculate Gravitation load for parallel robot
% P3PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G1P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:51
% EndTime: 2020-03-09 21:14:51
% DurationCPUTime: 0.24s
% Computational Cost: add. (387->70), mult. (359->117), div. (30->5), fcn. (300->18), ass. (0->55)
t509 = pkin(7) + qJ(2,3);
t499 = qJ(3,3) + t509;
t487 = sin(t499);
t490 = cos(t499);
t493 = sin(t509);
t496 = cos(t509);
t532 = 0.1e1 / (t487 * t496 - t493 * t490);
t510 = pkin(7) + qJ(2,2);
t500 = qJ(3,2) + t510;
t488 = sin(t500);
t491 = cos(t500);
t494 = sin(t510);
t497 = cos(t510);
t531 = 0.1e1 / (t488 * t497 - t494 * t491);
t511 = pkin(7) + qJ(2,1);
t501 = qJ(3,1) + t511;
t489 = sin(t501);
t492 = cos(t501);
t495 = sin(t511);
t498 = cos(t511);
t530 = 0.1e1 / (t489 * t498 - t495 * t492);
t516 = 0.1e1 / pkin(2);
t529 = 0.1e1 / pkin(3) * t516;
t512 = legFrame(3,3);
t502 = sin(t512);
t505 = cos(t512);
t481 = -t502 * g(1) + t505 * g(2);
t484 = t505 * g(1) + t502 * g(2);
t463 = -t490 * (mrSges(3,1) * t481 - mrSges(3,2) * t484) + (mrSges(3,1) * t484 + mrSges(3,2) * t481) * t487;
t508 = m(3) * pkin(2) + mrSges(2,1);
t528 = ((mrSges(2,2) * t484 - t481 * t508) * t496 + t493 * (t481 * mrSges(2,2) + t508 * t484) + t463) * t532 * t516;
t513 = legFrame(2,3);
t503 = sin(t513);
t506 = cos(t513);
t482 = -t503 * g(1) + t506 * g(2);
t485 = t506 * g(1) + t503 * g(2);
t464 = -t491 * (mrSges(3,1) * t482 - mrSges(3,2) * t485) + (mrSges(3,1) * t485 + mrSges(3,2) * t482) * t488;
t527 = ((mrSges(2,2) * t485 - t482 * t508) * t497 + t494 * (t482 * mrSges(2,2) + t508 * t485) + t464) * t531 * t516;
t514 = legFrame(1,3);
t504 = sin(t514);
t507 = cos(t514);
t483 = -t504 * g(1) + t507 * g(2);
t486 = t507 * g(1) + t504 * g(2);
t465 = -t492 * (mrSges(3,1) * t483 - mrSges(3,2) * t486) + (mrSges(3,1) * t486 + mrSges(3,2) * t483) * t489;
t526 = ((mrSges(2,2) * t486 - t483 * t508) * t498 + t495 * (t483 * mrSges(2,2) + t508 * t486) + t465) * t530 * t516;
t525 = t463 * t532 * t529;
t524 = t464 * t531 * t529;
t523 = t465 * t530 * t529;
t522 = t505 * t487 + t502 * t490;
t521 = -t487 * t502 + t505 * t490;
t520 = t506 * t488 + t503 * t491;
t519 = -t488 * t503 + t506 * t491;
t518 = t507 * t489 + t504 * t492;
t517 = -t489 * t504 + t507 * t492;
t1 = [t517 * t526 - (-pkin(2) * (t495 * t504 - t507 * t498) + t517 * pkin(3)) * t523 + t519 * t527 - (-pkin(2) * (t494 * t503 - t506 * t497) + t519 * pkin(3)) * t524 + t521 * t528 - (-pkin(2) * (t493 * t502 - t505 * t496) + t521 * pkin(3)) * t525 - g(1) * m(4); t518 * t526 - (pkin(2) * (t495 * t507 + t504 * t498) + t518 * pkin(3)) * t523 + t520 * t527 - (pkin(2) * (t494 * t506 + t503 * t497) + t520 * pkin(3)) * t524 + t522 * t528 - (pkin(2) * (t493 * t505 + t502 * t496) + t522 * pkin(3)) * t525 - g(2) * m(4); -(3 * g(3) * (m(1) + m(2) + m(3))) - g(3) * m(4);];
taugX  = t1;
