% Calculate Gravitation load for parallel robot
% P3PRRRR1G1P3A0
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:07
% EndTime: 2020-03-09 20:34:07
% DurationCPUTime: 0.25s
% Computational Cost: add. (177->66), mult. (323->130), div. (42->10), fcn. (264->18), ass. (0->60)
t585 = (m(1) + m(2) + m(3));
t633 = g(3) * t585;
t595 = legFrame(3,3);
t578 = sin(t595);
t581 = cos(t595);
t572 = -g(1) * t578 + g(2) * t581;
t605 = cos(qJ(3,3));
t589 = 0.1e1 / t605;
t599 = sin(qJ(3,3));
t575 = g(1) * t581 + g(2) * t578;
t600 = sin(qJ(2,3));
t606 = cos(qJ(2,3));
t617 = g(3) * t600 + t575 * t606;
t632 = ((mrSges(3,1) * t572 + t617 * mrSges(3,2)) * t605 + t599 * (t617 * mrSges(3,1) - mrSges(3,2) * t572)) * t589;
t596 = legFrame(2,3);
t579 = sin(t596);
t582 = cos(t596);
t573 = -g(1) * t579 + g(2) * t582;
t607 = cos(qJ(3,2));
t591 = 0.1e1 / t607;
t601 = sin(qJ(3,2));
t576 = g(1) * t582 + g(2) * t579;
t602 = sin(qJ(2,2));
t608 = cos(qJ(2,2));
t616 = g(3) * t602 + t576 * t608;
t631 = ((mrSges(3,1) * t573 + t616 * mrSges(3,2)) * t607 + t601 * (t616 * mrSges(3,1) - mrSges(3,2) * t573)) * t591;
t597 = legFrame(1,3);
t580 = sin(t597);
t583 = cos(t597);
t574 = -g(1) * t580 + g(2) * t583;
t609 = cos(qJ(3,1));
t593 = 0.1e1 / t609;
t603 = sin(qJ(3,1));
t577 = g(1) * t583 + g(2) * t580;
t604 = sin(qJ(2,1));
t610 = cos(qJ(2,1));
t615 = g(3) * t604 + t577 * t610;
t630 = ((mrSges(3,1) * t574 + t615 * mrSges(3,2)) * t609 + t603 * (t615 * mrSges(3,1) - mrSges(3,2) * t574)) * t593;
t586 = 0.1e1 / t600;
t629 = t586 * t589;
t587 = 0.1e1 / t602;
t628 = t587 * t591;
t588 = 0.1e1 / t604;
t627 = t588 * t593;
t626 = t599 * t606;
t625 = t601 * t608;
t624 = t603 * t610;
t623 = t605 * t606;
t622 = t607 * t608;
t621 = t609 * t610;
t598 = mrSges(2,2) - mrSges(3,3);
t584 = g(3) * t598;
t614 = mrSges(3,1) * t605 - mrSges(3,2) * t599 + mrSges(2,1);
t620 = (t600 * t584 + (t598 * t606 + t600 * t614) * t575 - t614 * t606 * g(3)) * t586 / t605 ^ 2;
t613 = mrSges(3,1) * t607 - mrSges(3,2) * t601 + mrSges(2,1);
t619 = (t602 * t584 + (t598 * t608 + t602 * t613) * t576 - t613 * t608 * g(3)) * t587 / t607 ^ 2;
t612 = mrSges(3,1) * t609 - mrSges(3,2) * t603 + mrSges(2,1);
t618 = (t604 * t584 + (t598 * t610 + t604 * t612) * t577 - t612 * t610 * g(3)) * t588 / t609 ^ 2;
t611 = 0.1e1 / pkin(2);
t1 = [-g(1) * m(4) + (-(t580 * t603 + t583 * t621) * t627 - (t579 * t601 + t582 * t622) * t628 - (t578 * t599 + t581 * t623) * t629) * t633 + ((-t580 * t624 - t583 * t609) * t618 + t580 * t630 + (-t579 * t625 - t582 * t607) * t619 + t579 * t631 + (-t578 * t626 - t581 * t605) * t620 + t578 * t632) * t611; -g(2) * m(4) + (-(t580 * t621 - t583 * t603) * t627 - (t579 * t622 - t582 * t601) * t628 - (t578 * t623 - t581 * t599) * t629) * t633 + ((-t580 * t609 + t583 * t624) * t618 - t583 * t630 + (-t579 * t607 + t582 * t625) * t619 - t582 * t631 + (-t578 * t605 + t581 * t626) * t620 - t581 * t632) * t611; (-m(4) - (3 * t585)) * g(3);];
taugX  = t1;
